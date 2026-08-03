#' Reconcile an estimated ie with the numerical floor used to build the covariance matrix
#'
#' spcov_matrix.*() floors ie when building the covariance matrix and reports the value back.
#' Only estimated ie are
#' touched -- known/fixed ie (including spcov_type = "none") is left alone.
#'
#' @noRd
floor_estimated_ie <- function(spcov_params_val, spcov_is_known, diagtol) {
  # car/sar (areal) covariance never floors ie -- spcov_matrix.car()/.sar()
  # ignore diagtol entirely -- so there is nothing to reconcile
  if (inherits(spcov_params_val, c("car", "sar"))) {
    return(spcov_params_val)
  }
  # only reconcile ie that was actually estimated by the optimizer; a known/fixed
  # ie (spcov_type = "none"'s type-mandated ie = 0, or a user-supplied known value
  # via spcov_initial()) is left exactly as reported
  if (isTRUE(spcov_is_known[["ie"]])) {
    return(spcov_params_val)
  }
  de <- if ("de" %in% names(spcov_params_val)) spcov_params_val[["de"]] else 0
  # matches the exact floor spcov_matrix.*() applies when building the covariance matrix
  spcov_params_val[["ie"]] <- max(spcov_params_val[["ie"]], 1e-4 * de, diagtol)
  spcov_params_val
}
