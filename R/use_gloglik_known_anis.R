# use_gloglik* family overview: see use_gloglik_known.R for the two axes that
# distinguish these sibling files (anisotropy, known-vs-estimated parameters).
# This file: anisotropic + known parameters. Like use_gloglik_known.R, no
# optim() search happens -- but unlike it, the distance matrix must still be
# built here (rather than passed in) because it depends on the fixed
# rotate/scale anisotropy parameters.
#' Use Gaussian log-likelihood estimation with anisotropy and known covariance parameters
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param estmethod The estimation method (\code{"reml"} or \code{"ml"})
#' @param X Model matrix
#' @param y Response vector
#' @param n Sample size
#' @param p Number of fixed effects
#' @param xcoord_val x-coordinate value
#' @param ycoord_val y-coordinate value
#' @param randcov_initial A \code{randcov_initial} object
#' @param randcov_Zs Random effects design matrices
#' @param observed_index The index of observed values
#' @param partition_matrix The partition matrix
#'
#' @return Known covariance parameters
#'
#' @noRd
use_gloglik_known_anis <- function(spcov_initial, data_object, estmethod, randcov_initial) {
  spcov_params_val <- get_spcov_params(class(spcov_initial), spcov_initial$initial)
  randcov_params_val <- randcov_params(randcov_initial$initial)

  # rotate/rescale coordinates by the fixed anisotropy parameters before
  # computing distances, since anisotropic covariance is isotropic in this
  # transformed coordinate space
  dist_matrix_list <- build_anis_dist_matrix_list(spcov_params_val, data_object)
  ## find relevant products
  gll_prods <- gloglik_products(
    spcov_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )
  ## compute -2ll
  minustwologlik <- get_minustwologlik(gll_prods, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE)

  # return parameter values and optim output
  optim_output <- known_optim_output_stub(minustwologlik)

  # return list
  list(
    spcov_params_val = spcov_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, randcov = randcov_initial$is_known)
  )
}
