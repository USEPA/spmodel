get_vcov_glm <- function(cov_betahat) {
  vcov_fixed <- cov_betahat
  if (any(diag(vcov_fixed) < 0)) {
    warning("Model fit potentially unstable. Consider fixing ie (via spcov_initial) at some non-zero value greater than 1e-4 and refitting the model.", call. = FALSE)
  }
  vcov <- list(fixed = vcov_fixed)
}
