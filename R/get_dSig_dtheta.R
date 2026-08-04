#' @noRd
dSig_dtheta_spcov <- function(spcov_params_val, dist_matrix, ...) {
  UseMethod("dSig_dtheta_spcov", spcov_params_val)
}

#' @noRd
dSig_dtheta_spcov.default <- function(spcov_params_val, dist_matrix, ...) {
  stop(
    "Closed-form covariance parameter derivatives are not implemented for spcov_type = \"",
    class(spcov_params_val)[1], "\" yet; use method = \"numeric\" instead.",
    call. = FALSE
  )
}

#' @noRd
dSig_dtheta_spcov.exponential <- function(spcov_params_val, dist_matrix, ...) {
  de <- spcov_params_val[["de"]]
  range <- spcov_params_val[["range"]]
  R <- exp(-dist_matrix / range)
  list(
    de = R,
    ie = Diagonal(nrow(dist_matrix)),
    range = de * R * (dist_matrix / range^2)
  )
}

#' @noRd
dSig_dtheta_randcov <- function(randcov_Zs, randcov_names) {
  val <- lapply(randcov_names, function(nm) as.matrix(randcov_Zs[[nm]][["ZZt"]]))
  names(val) <- randcov_names
  val
}

get_dSig_dtheta_cov <- function(context) {

  dSig_spcov <- dSig_dtheta_spcov(context$spcov_params, context$dist_matrix)[context$spcov_names_free]
  dSig_randcov <- if (!is.null(context$randcov_names_free)) {
    dSig_dtheta_randcov(context$randcov_Zs, context$randcov_names_free)
  } else {
    list()
  }
  c(dSig_spcov, dSig_randcov)
}
