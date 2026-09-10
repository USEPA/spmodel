#' Evaluate minus twice the Laplace-approximated log-likelihood for \code{optim()}
#'
#' @param par The current optimization parameter vector
#' @param spcov_orig2optim A spatial covariance parameter list on the optimization scale
#' @param dispersion_orig2optim A dispersion parameter list on the optimization scale
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param dist_matrix_list A list of distance matrices
#' @param spcov_profiled Whether the spatial covariance parameters are profiled out
#' @param randcov_orig2optim A random effect variance list on the optimization
#'   scale (or \code{NULL} if there are no random effects)
#' @param randcov_profiled Whether the random effect variances are profiled out
#'
#' @return Minus twice the Laplace-approximated (restricted) log-likelihood at
#'   \code{par}, after transforming \code{par} back to the original covariance
#'   parameter scales
#'
#' @noRd
laploglik <- function(par, spcov_orig2optim, dispersion_orig2optim, data_object, estmethod, dist_matrix_list,
                      spcov_profiled, randcov_orig2optim = NULL,
                      randcov_profiled = NULL) {
  # optim() passes a single flat parameter vector; the dispersion element(s)
  # are pulled off the front first and converted back to their original
  # scale, leaving the remaining (spatial/random effect) parameters in par
  # dispersion first then remove
  dispersion_result <- peel_dispersion(dispersion_orig2optim, par, data_object$family)
  dispersion_params_val <- dispersion_result$dispersion_params_val
  par <- dispersion_result$par

  unpacked <- unpack_optim2orig(spcov_orig2optim, randcov_orig2optim, par, spcov_profiled, randcov_profiled, data_object)
  spcov_params_val <- unpacked$spcov_params_val
  randcov_params_val <- unpacked$randcov_params_val

  lapll_prods <- laploglik_products(
    spcov_params_val, dispersion_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )

  # spcov_profiled/randcov_profiled are forced to FALSE here because, unlike
  # the Gaussian case, the GLM Laplace approximation re-optimizes the latent
  # w vector each call, so profiling out variances in closed form does not apply
  minustwolaploglik <- get_minustwolaploglik(lapll_prods, estmethod, data_object$n,
    data_object$p,
    spcov_profiled = FALSE, randcov_profiled = FALSE
  )
}
