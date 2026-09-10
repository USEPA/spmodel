#' Get random effects design matrices for every random effect term
#'
#' @param data The data
#' @param randcov_names The names of the random effects (as returned by \code{get_randcov_names()})
#' @param ZZt Should \eqn{ZZ'} be returned for each term?
#' @param ZtZ Should \eqn{Z'Z} be returned for each term?
#' @param xlev_list Optional list of factor levels to enforce for each random
#'   effect term (used so a design matrix built for a data subset, e.g. a
#'   big-data local neighborhood, still spans the same levels as the full data)
#'
#' @return A named list (one element per random effect term) of design matrices
#'
#' @noRd
get_randcov_Zs <- function(data, randcov_names = NULL, ZZt = TRUE, ZtZ = FALSE, xlev_list = NULL) {
  if (is.null(randcov_names)) {
    randcov_Zs <- NULL
  } else {
    randcov_Zs <- lapply(randcov_names, get_randcov_Z, data, ZZt, ZtZ, xlev_list)
    names(randcov_Zs) <- randcov_names
  }
  randcov_Zs
}

#' Get the random effects design matrix for a single random effect term
#'
#' @param randcov_name The name of a single random effect term (as returned by \code{get_randcov_names()})
#' @param data The data
#' @param ZZt Should \eqn{ZZ'} be returned?
#' @param ZtZ Should \eqn{Z'Z} be returned?
#' @param xlev_list Optional list of factor levels to enforce for each random
#'   effect term (used so a design matrix built for a data subset, e.g. a
#'   big-data local neighborhood, still spans the same levels as the full data)
#'
#' @return A list with elements \code{Z}, \code{ZZt}, and \code{ZtZ}
#'
#' @noRd
get_randcov_Z <- function(randcov_name, data, ZZt = TRUE, ZtZ = FALSE, xlev_list = NULL) {
  bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
  Z_reform <- reformulate(bar_split[[2]], intercept = FALSE)
  if (is.null(xlev_list)) {
    Z_frame <- model.frame(Z_reform, data = data, drop.unused.levels = FALSE)
  } else {
    Z_frame <- model.frame(Z_reform, data = data, drop.unused.levels = FALSE, xlev = xlev_list[[randcov_name]])
  }
  if (any(!attr(terms(Z_frame), "dataClasses") %in% c("character", "factor", "ordered"))) {
    stop("Random effect grouping variables must be categorical or factor.", call. = FALSE)
  }
  # built directly as a sparse indicator from factor codes rather than via a
  # dense model.matrix(), which would allocate an n x nlevels dense matrix
  # before ever going sparse -- expensive when the grouping variable has many levels
  Z_factor <- Z_frame[[1]]
  if (!is.factor(Z_factor)) Z_factor <- factor(Z_factor)
  Z_index <- Matrix::t(Matrix::fac2sparse(Z_factor, drop.unused.levels = FALSE))
  # fac2sparse() names columns by bare factor level (e.g. "a"), but downstream
  # per-level output (e.g. fitted(object, type = "randcov")) and the
  # model.matrix()-based prediction-time equivalent (get_randcov_vectors())
  # both expect model.matrix()'s "varname + level" convention (e.g. "groupa")
  # -- restore it here so level names match between fitting and prediction
  colnames(Z_index) <- get_factor_level_names(names(Z_frame)[1], levels(Z_factor))
  if (bar_split[[1]] == "1") {
    # "1 | group" (random intercept): Z is just the grouping indicator matrix
    Z <- Z_index
  } else {
    # "x | group" (random slope): each row of the grouping indicator matrix
    # is additionally scaled by that row's x value, so the random effect for
    # a group applies to the slope on x rather than a flat intercept shift
    Z_mod_reform <- reformulate(bar_split[[1]], intercept = FALSE)
    # don't need xlev here as Z_mod_reform is for the continuous variable x
    Z_mod_frame <- model.frame(Z_mod_reform, data = data, drop.unused.levels = FALSE)
    Z_mod <- model.matrix(Z_mod_reform, Z_mod_frame)
    if (NCOL(Z_mod) > 1) {
      stop("All variable names to the left of | in random must be numeric.", call. = FALSE)
    }
    Z <- as.vector(Z_mod) * Z_index
  }
  if (ZZt) {
    ZZt <- tcrossprod(Z, Z)
  } else {
    ZZt <- NULL
  }

  if (ZtZ) {
    ZtZ <- crossprod(Z, Z)
  } else {
    ZtZ <- NULL
  }

  list(Z = Z, ZZt = ZZt, ZtZ = ZtZ)
}

#' Get each observation's random effect group label and slope coefficient,
#' one random effect term at a time, without ever forming \code{ZZt}
#'
#' Large data (efficient) companion to \code{get_randcov_Zs()}/\code{get_randcov_Z()}:
#' returns, per random effect term, the length-n group-label vector and
#' length-n "coefficient" vector (all \code{1} for a random intercept, the
#' numeric predictor's value for a random slope) instead of the full n x n
#' \code{ZZt} matrix. One random
#' effect term's (i, j) covariance contribution is
#' \code{variance * coef_i * coef_j * (group_i == group_j)}; see
#' \code{get_randcov_local()}, which assembles this (summed across terms)
#' for a small set of indices.
#'
#' @param random The random effect formula, or \code{NULL}.
#' @param data Data.
#'
#' @return A named list (one element per \code{get_randcov_names(random)}
#'   term) of lists with elements \code{group} (factor) and \code{coef}
#'   (numeric, length \code{NROW(data)}), or \code{NULL} if \code{random} is
#'   \code{NULL}.
#'
#' @noRd
get_randcov_groups <- function(random, data) {
  if (is.null(random)) {
    return(NULL)
  }
  randcov_names <- get_randcov_names(random)
  groups <- lapply(randcov_names, function(randcov_name) {
    bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
    Z_reform <- reformulate(bar_split[[2]], intercept = FALSE)
    Z_frame <- model.frame(Z_reform, data = data, drop.unused.levels = FALSE)
    if (any(!attr(terms(Z_frame), "dataClasses") %in% c("character", "factor", "ordered"))) {
      stop("Random effect grouping variables must be categorical or factor.", call. = FALSE)
    }
    group <- Z_frame[[1]]
    if (!is.factor(group)) group <- factor(group)
    if (bar_split[[1]] == "1") {
      # "1 | group" (random intercept): every row contributes coefficient 1
      coef <- rep(1, NROW(data))
    } else {
      # "x | group" (random slope): matches get_randcov_Z()'s Z_mod handling
      Z_mod_reform <- reformulate(bar_split[[1]], intercept = FALSE)
      Z_mod_frame <- model.frame(Z_mod_reform, data = data, drop.unused.levels = FALSE)
      Z_mod <- model.matrix(Z_mod_reform, Z_mod_frame)
      if (NCOL(Z_mod) > 1) {
        stop("All variable names to the left of | in random must be numeric.", call. = FALSE)
      }
      coef <- as.vector(Z_mod)
    }
    list(group = group, coef = coef)
  })
  names(groups) <- randcov_names
  groups
}

#' Build a small random effect covariance vector/matrix contribution from
#' \code{get_randcov_groups()}'s efficient per-term group/coefficient
#' vectors, for an arbitrary pair of index sets
#'
#' The (i, j) contribution of one random effect term is
#' \code{variance * coef_i * coef_j * (group_i == group_j)}; contributions
#' from every term (named in \code{randcov_params}) are summed, matching
#' \code{randcov_matrix()}'s \code{Reduce("+", ...)} across terms.
#'
#' @param randcov_params A named numeric vector of random effect variances
#'   (names matching \code{get_randcov_names(random)}), or \code{NULL}.
#' @param randcov_groups The named list from \code{get_randcov_groups()}, or
#'   \code{NULL}.
#' @param idx1,idx2 Row indices into each term's \code{group}/\code{coef}
#'   vectors.
#'
#' @return A matrix of dimension \code{length(idx1) x length(idx2)}, or
#'   \code{NULL} if either input is \code{NULL}.
#'
#' @noRd
get_randcov_local <- function(randcov_params, randcov_groups, idx1, idx2) {
  if (is.null(randcov_params) || is.null(randcov_groups)) {
    return(NULL)
  }
  contribs <- lapply(names(randcov_params), function(term) {
    g <- randcov_groups[[term]]
    same_group <- outer(g$group[idx1], g$group[idx2], FUN = "==")
    coef_prod <- outer(g$coef[idx1], g$coef[idx2], FUN = "*")
    randcov_params[[term]] * same_group * coef_prod
  })
  Reduce("+", contribs)
}

#' Get the factor levels of a single random effect term's grouping variable
#'
#' @param randcov_name The name of a single random effect term (as returned by \code{get_randcov_names()})
#' @param data The data
#'
#' @return The factor levels of the grouping variable, as returned by \code{.getXlevels()}
#'
#' @noRd
get_randcov_xlev <- function(randcov_name, data) {
  # used to capture the full-data factor levels up front, so they can later
  # be passed as xlev_list into get_randcov_Z()/get_randcov_Zs() and enforced
  # on a data subset (e.g. a local-estimation neighborhood or newdata at
  # prediction time) that might not itself contain every level
  bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
  Z_reform <- reformulate(bar_split[[2]], intercept = FALSE)
  Z_frame <- model.frame(Z_reform, data = data, drop.unused.levels = FALSE)
  .getXlevels(terms(Z_frame), Z_frame)
}
