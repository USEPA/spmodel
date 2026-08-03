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
  # bundle the three coefficient types into one list so coef() can dispatch
  # on its type argument ("fixed", "spcov", "randcov") without recomputation
  list(fixed = betahat, spcov = spcov_params, randcov = randcov_params)
}
