#' Get variance covariance matrices
#'
#' @param cov_betahat The covariance of the fixed effects
#'
#' @return Variance covariance matrices
#'
#' @noRd
get_vcov <- function(cov_betahat) {
  ## variance covariance matrix
  vcov_fixed <- cov_betahat
  vcov_spcov <- NULL
  vcov_randcov <- NULL
  vcov <- list(fixed = vcov_fixed, spcov = vcov_spcov, randcov = vcov_randcov)
}
