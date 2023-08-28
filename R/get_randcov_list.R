#' Get random effects design matrices
#'
#' @param randcov_vars The names of the random effects
#' @param data The data
#' @param ZZt Should ZZt be returned?
#' @param ZtZ Should ZtZ be returned?
#'
#' @return Random effects design matrices
#'
#' @noRd
get_randcov_list <- function(index, randcov_Zs = NULL, randcov_names = NULL) {
  if (is.null(randcov_names)) {
    randcov_list <- NULL
  } else {
    index_vals <- split(seq_along(index), index)
    randcov_list <- lapply(index_vals, function(x) get_randcov_index_list(x, index, randcov_Zs, randcov_names))
    names(randcov_list) <- names(index_vals)
  }
  randcov_list
}

get_randcov_index_list <- function(index_val, index, randcov_Zs, randcov_names) {
  Z_lists <- lapply(randcov_names, function(x) get_randcov_var_list(x, index_val, index, randcov_Zs))
  names(Z_lists) <- randcov_names
  Z_lists
}

get_randcov_var_list <- function(randcov_name, index_val, index, randcov_Zs) {
  Z_list <- randcov_Zs[[randcov_name]][["Z"]][index_val, , drop = FALSE]
  if (is.null(randcov_Zs[[randcov_name]][["ZZt"]])) {
    ZZt_list <- NULL
  } else {
    ZZt_list <- randcov_Zs[[randcov_name]][["ZZt"]][index_val, index_val, drop = FALSE]
  }
  if (is.null(randcov_Zs[[randcov_name]][["ZtZ"]])) {
    ZtZ_list <- NULL
  } else {
    ZtZ_list <- randcov_Zs[[randcov_name]][["ZtZ"]][index_val, index_val, drop = FALSE]
  }
  list(Z = Z_list, ZZt = ZZt_list, ZtZ = ZtZ_list)
}

get_randcov_Zs <- function(data, randcov_names = NULL, ZZt = TRUE, ZtZ = FALSE, xlev_list = NULL) {
  if (is.null(randcov_names)) {
    randcov_Zs <- NULL
  } else {
    randcov_Zs <- lapply(randcov_names, get_randcov_Z, data, ZZt, ZtZ, xlev_list)
    names(randcov_Zs) <- randcov_names
  }
  randcov_Zs
}

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
  Z_index <- Matrix(model.matrix(Z_reform, Z_frame), sparse = TRUE)
  if (bar_split[[1]] == "1") {
    Z <- Z_index
  } else {
    Z_mod_reform <- reformulate(bar_split[[1]], intercept = FALSE)
    Z_mod_frame <- model.frame(Z_mod_reform, data = data, drop.unused.levels = FALSE)
    Z_mod <- model.matrix(Z_mod_reform, Z_mod_frame)
    if (NCOL(Z_mod) > 1) {
      stop("All variable names to the left of | in random must be numeric.", call. = FALSE)
    }
    Z <- as.vector(Z_mod) * Z_index
  }
  # find and replace values not observed
  Z_levels_observed <- which(colSums(abs(Z)) > 0)
  Z <- Z[, Z_levels_observed, drop = FALSE]
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
