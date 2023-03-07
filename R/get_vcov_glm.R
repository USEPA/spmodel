get_vcov_glm <- function(cov_betahat) {
  vcov_fixed <- cov_betahat
  vcov <- list(fixed = vcov_fixed)
}
