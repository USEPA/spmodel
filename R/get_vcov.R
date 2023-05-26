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
  if (any(diag(vcov_fixed) < 0)) {
    warning("Model fit potentially unstable. Consider fixing ie (via spcov_initial) at some non-zero value greater than 1e-4 and refitting the model.", call. = FALSE)
  }
  vcov_spcov <- NULL
  vcov_randcov <- NULL
  vcov <- list(fixed = vcov_fixed, spcov = vcov_spcov, randcov = vcov_randcov)
}
