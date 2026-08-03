#' Perform initial checks for spautor models
#'
#' @param spcov_type The spatial covariance type
#' @param y_coord_given Is the y-coordinate given?
#'
#' @return Error messages or nothing
#'
#' @noRd
spglm_checks <- function(family, spcov_initial, xcoord_given, ycoord_given, estmethod, anisotropy, random_given) {
  spcov_type <- class(spcov_initial)
  check_not_areal_type(spcov_type, "spglm", "spgautor")
  warn_ycoord_ignored_for_1d_cov(spcov_type, ycoord_given)
  check_estmethod_reml_ml(estmethod, "Estimation method must be \"reml\" or \"ml\".")
  check_family_valid(family)
}
