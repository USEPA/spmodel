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
  # optimizers work on an unconstrained/transformed ("optim") parameter scale
  # (e.g. log-transformed variances), so map par back to the original,
  # interpretable covariance parameter scale before evaluating the likelihood
  unpacked <- unpack_optim2orig(spcov_orig2optim, randcov_orig2optim, par, spcov_profiled, randcov_profiled, data_object)
  spcov_params_val <- unpacked$spcov_params_val
  randcov_params_val <- unpacked$randcov_params_val

  gll_prods <- gloglik_products(
    spcov_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )

  # combines the covariance-parameter likelihood products into the final
  # -2*log-likelihood (or REML/other objective) value the optimizer minimizes
  minustwologlik <- get_minustwologlik(gll_prods, estmethod, data_object$n, data_object$p, spcov_profiled, randcov_profiled = randcov_profiled)
}
