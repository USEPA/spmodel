#' Get the minus twice log likelihood
#'
#' @param gloglik_products A \code{gloglik_products} object
#' @param estmethod The estimation method (ml or reml)
#' @param n The sample size
#' @param p The number of fixed effects
#' @param spcov_profiled Is the overall variance of the spatial covariance profiled?
#'
#' @return minus twice the log likelihood
#'
#' @noRd
get_minustwologlik <- function(gloglik_products, estmethod, n, p,
                               spcov_profiled, randcov_profiled = NULL) {
  # when the overall variance (sigma^2) is profiled out, it has an analytic
  # solution given the other covariance parameters, so the likelihood formula
  # substitutes that solution in directly (log(l2) terms below) instead of
  # optimizing over sigma^2 as a free parameter -- this shrinks the
  # optimization problem by one dimension
  if (spcov_profiled && (is.null(randcov_profiled) || randcov_profiled)) {
    if (estmethod == "reml") {
      minustwologlik <- as.numeric(gloglik_products$l1 + (n - p) * log(gloglik_products$l2) + gloglik_products$l3 + (n - p) * (1 + log(2 * pi / (n - p))))
    } else if (estmethod == "ml") {
      minustwologlik <- as.numeric(gloglik_products$l1 + n * log(gloglik_products$l2) + n * (1 + log(2 * pi / n)))
    }
  } else {
    # sigma^2 not profiled out -- l2 already incorporates it directly
    if (estmethod == "reml") {
      minustwologlik <- as.numeric(gloglik_products$l1 + gloglik_products$l2 + gloglik_products$l3 + (n - p) * log(2 * pi))
    } else if (estmethod == "ml") {
      minustwologlik <- as.numeric(gloglik_products$l1 + gloglik_products$l2 + n * log(2 * pi))
    }
  }
}
