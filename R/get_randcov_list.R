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
