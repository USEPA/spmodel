#' Use Gaussian log-likelihood estimation when covariance parameters are known
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param estmethod The estimation method (\code{"reml"} or \code{"ml"})
#' @param X Model matrix
#' @param y Response vector
#' @param n Sample size
#' @param p Number of fixed effects
#' @param dist_matrix Distance matrix (Euclidean or neighbors)
#' @param randcov_initial A \code{randcov_initial} object
#' @param randcov_Zs Random effects design matrices
#' @param observed_index The index of observed values
#' @param partition_matrix The partition matrix
#'
#' @return Known covariance parameters
#'
#' @noRd
use_gloglik_known <- function(spcov_initial, data_object, estmethod, dist_matrix_list, randcov_initial) {
  spcov_params_val <- get_spcov_params(class(spcov_initial), spcov_initial$initial)

  randcov_params_val <- randcov_params(randcov_initial$initial)
  ## find relevant products
  gll_prods <- gloglik_products(
    spcov_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )
  ## compute -2ll
  minustwologlik <- get_minustwologlik(gll_prods, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE)
  # return parameter values and optim output
  optim_output <- list(
    method = NA, control = NA, value = minustwologlik,
    counts = NA, convergence = NA,
    message = NA, hessian = NA
  )

  # return list
  list(
    spcov_params_val = spcov_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, randcov = randcov_initial$is_known)
  )
}
