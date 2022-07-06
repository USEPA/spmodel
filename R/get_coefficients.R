#' Organize relevant coefficients
#'
#' @param betahat Fixed effects
#' @param spcov_params Spatial covariance parameters
#' @param randcov_params Random effects
#'
#' @return A list of relevant coefficients
#'
#' @noRd
get_coefficients <- function(betahat, spcov_params, randcov_params = NULL) {
  list(fixed = betahat, spcov = spcov_params, randcov = randcov_params)
}
