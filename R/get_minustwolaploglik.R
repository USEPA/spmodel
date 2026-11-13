#' Compute minus twice the Laplace-approximated log-likelihood
#'
#' @param laploglik_products A \code{laploglik_products} object
#' @param estmethod The estimation method (\code{"reml"} or \code{"ml"})
#' @param n The sample size
#' @param p The number of fixed effects
#' @param spcov_profiled Whether the spatial covariance parameters are profiled out (currently always \code{FALSE})
#' @param randcov_profiled Whether the random effect variances are profiled out (currently always \code{FALSE})
#'
#' @return Minus twice the Laplace-approximated (restricted) log-likelihood,
#'   the objective minimized during covariance parameter estimation for
#'   \code{spglm()} and \code{spgautor()} models
#'
#' @noRd
get_minustwolaploglik <- function(laploglik_products, estmethod, n, p, spcov_profiled = FALSE, randcov_profiled = FALSE) {
# On the -2log-likelihood scale, the Laplace constant cancels the Gaussian n log(2 pi)
  if (estmethod == "reml") {
    minustwolaploglik <- as.numeric(laploglik_products$l00 + laploglik_products$l01 + laploglik_products$l1 +
      laploglik_products$l2 + laploglik_products$l3 +
      -p * log(2 * pi))
  } else if (estmethod == "ml") {
    minustwolaploglik <- as.numeric(laploglik_products$l00 + laploglik_products$l01 + laploglik_products$l1 +
      laploglik_products$l2)
  }

  minustwolaploglik
}
