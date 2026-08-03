# use_gloglik* family overview: these functions evaluate/optimize the Gaussian
# log-likelihood for spatial linear models (splm/spautor). Variants differ along
# two axes -- anisotropy (plain vs "_anis") and whether covariance parameters are
# already known (plain vs "_known"). This file: isotropic + known parameters, so
# no optim() search is run -- the fixed parameters are simply plugged in and the
# resulting -2*loglik is returned. Compare with use_gloglik.R (isotropic, estimated),
# use_gloglik_anis.R (anisotropic, estimated), and use_gloglik_known_anis.R
# (anisotropic, known).
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
  # optim() was never called (parameters are fixed), so its diagnostic fields
  # are filled with NA -- only value (the -2ll) is meaningful here
  optim_output <- known_optim_output_stub(minustwologlik)

  # return list
  list(
    spcov_params_val = spcov_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, randcov = randcov_initial$is_known)
  )
}
