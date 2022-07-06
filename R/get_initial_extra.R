#' Get initial values for the extra parameter (in four parameter geostatistical models)
#'
#' @param spcov_type The spatial covariance type
#'
#' @return An extra initial value
#'
#' @noRd
get_initial_extra <- function(spcov_type) {
  if (spcov_type == "matern") {
    initial_extra <- 2
  } else if (spcov_type == "cauchy") {
    initial_extra <- 1
  } else if (spcov_type == "pexponential") {
    initial_extra <- 0.8
  } else {
    initial_extra <- NULL # so it can be run with 3 parameter models without a logical statement
  }
  initial_extra
}
