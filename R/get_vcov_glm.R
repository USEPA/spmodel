get_vcov_glm <- function(cov_betahat_corrected, cov_betahat_uncorrected) {
  if (any(diag(cov_betahat_corrected) < 0)) {
    warning("Model fit potentially unstable. Consider fixing ie (via spcov_initial) at some non-zero value greater than 1e-4 and refitting the model.", call. = FALSE)
  }
  vcov_fixed <- list(corrected = cov_betahat_corrected, uncorrected = cov_betahat_uncorrected)
  vcov <- list(fixed = vcov_fixed)
}
