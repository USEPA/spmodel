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
  # this needs some work for now -- three reasons
  # 1. it needs to not include unknown parameters in the inverse covariance parameter covariance matrix
  # 2. it needs to use the original covariance parameters when estimation was done using profiling
  # 3. it needs to use the non-optimized scale (or mann-wald from the optimized to the non-optimized)
  # if (dotlist$optim$hessian && estmethod %in% c("reml", "ml")) {
  #   # hess is minus fisher information
  #   # min -2 log lik gives partials matrix of minus 2 l(theta) (hessian)
  #   # divide by 2 gives partial of minus l(theta) which is information
  #   # solve information (minus hessian / 2)
  #   halfhessian_val <- cov_params_estimated$optim_output$hessian / 2
  #   diag(halfhessian_val) <- diag(halfhessian_val) + 1e-4
  #   vcov_covariance <- solve(halfhessian_val)
  # } else {
  #   vcov_covariance <- NULL
  # }
  vcov <- list(fixed = vcov_fixed, spcov = vcov_spcov, randcov = vcov_randcov)
}
