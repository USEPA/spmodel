#' Find (minus twice negative) Gaussian log-likelihood while optimizing
#'
#' @param par Parameters to optimize over
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param estmethod The estimation method
#' @param X Model matrix
#' @param y Response vector
#' @param n Sample size
#' @param p Number of fixed effects
#' @param dist_matrix Distance matrix (Euclidean or neighbor)
#' @param spcov_profiled Whether the overall spatial variance is profiled
#' @param randcov_orig2optim A \code{randcov_orig2optim} object
#' @param randcov_Zs Random effect design matrices
#' @param observed_index Index of observed values
#' @param partition_matrix Partition matrix
#'
#' @return (Minus twice negative) Gaussian log-likelihood
#'
#' @noRd
gloglik <- function(par, spcov_orig2optim, data_object, estmethod, dist_matrix_list,
                    spcov_profiled, randcov_orig2optim = NULL,
                    randcov_profiled = NULL) {


  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim, par, spcov_profiled = spcov_profiled, data_object = data_object)

  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim), spcov_orig_val = spcov_orig_val)

  #
  # transforming to original scale
  randcov_orig_val <- randcov_optim2orig(randcov_orig2optim, spcov_orig2optim, par,
    randcov_profiled = randcov_profiled,
    spcov_optim2orig = spcov_params_val
  )

  # need to deal with list if randcov_profiled as sp variance changes
  if (!is.null(randcov_profiled) && randcov_profiled) {
    spcov_params_val <- randcov_orig_val$spcov_optim2orig
    randcov_orig_val <- randcov_orig_val$fill_orig_val
  }

  # making a random effects vector
  randcov_params_val <- randcov_params(randcov_orig_val)
  gll_prods <- gloglik_products(
    spcov_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )

  minustwologlik <- get_minustwologlik(gll_prods, estmethod, data_object$n, data_object$p, spcov_profiled, randcov_profiled = randcov_profiled)
}
