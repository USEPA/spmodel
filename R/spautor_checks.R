#' Perform initial checks for spautor models
#'
#' @param spcov_type The spatial covariance type
#' @param W_given Is the spatial weight matrix given?
#' @param data data
#'
#' @return Error messages or nothing
#'
#' @noRd
spautor_checks <- function(spcov_type, W_given, data, estmethod) {
  check_not_point_referenced_type(spcov_type)
  check_W_given_data_class(W_given, data)
  check_estmethod_reml_ml(estmethod, "Estimation method must be \"reml\", or \"ml\".")
}
