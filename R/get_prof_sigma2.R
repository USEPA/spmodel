#' Get overall variance if profiling used
#'
#' @param spcov_params_val A \code{spcov_params} object
#' @param dist_matrix A distance matrix
#' @param estmethod The estimation method
#' @param X A model matrix
#' @param y A response vector
#' @param n The sample size
#' @param p The number of fixed effects
#' @param observed_index A vector indicating observed values
#'
#' @return The overall variance
#'
#' @noRd
get_prof_sigma2 <- function(spcov_params_val, data_object, estmethod,
                            dist_matrix_list, randcov_params_val = NULL) {
  gll_prods <- gloglik_products(
    spcov_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )
  ## finding sigma2 estimate (as we are have a spcov_profiled likelihood)
  if (estmethod == "reml") {
    sigma2 <- gll_prods$l2 / (data_object$n - data_object$p)
  } else if (estmethod == "ml") {
    sigma2 <- gll_prods$l2 / data_object$n
  }
  as.numeric(sigma2)
}
