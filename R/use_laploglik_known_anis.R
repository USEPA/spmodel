# see use_laploglik.R for an overview of the use_laploglik* family; this file
# is the anisotropic + known-parameters variant (no optim() search, but the
# distance matrix must still be rebuilt from the fixed rotate/scale values).
#' Evaluate the Laplace-approximated log-likelihood at fully-known parameters under anisotropy
#'
#' @param spcov_initial A \code{spcov_initial} object (all parameters fixed)
#' @param dispersion_initial A \code{dispersion_initial} object (fixed)
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param dist_matrix_list A list of distance matrices (unused; the anisotropy-corrected
#'   distances are computed internally from the fixed rotate/scale parameters)
#' @param randcov_initial A \code{randcov_initial} object (fixed, or \code{NULL})
#'
#' @return The same value as \code{use_laploglik_known()}, using
#'   anisotropy-corrected distances
#'
#' @noRd
use_laploglik_known_anis <- function(spcov_initial, dispersion_initial, data_object, estmethod, dist_matrix_list, randcov_initial) {
  spcov_params_val <- get_spcov_params(class(spcov_initial), spcov_initial$initial)
  dispersion_params_val <- dispersion_params(data_object$family, dispersion_initial$initial)
  randcov_params_val <- randcov_params(randcov_initial$initial)

  # rotate/rescale coordinates by the fixed anisotropy parameters so the
  # resulting distances make the covariance isotropic in transformed space
  dist_matrix_list <- build_anis_dist_matrix_list(spcov_params_val, data_object)

  lapll_prods <- laploglik_products(
    spcov_params_val, dispersion_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )
  ## compute -2ll
  minustwolaploglik <- get_minustwolaploglik(lapll_prods, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE)
  # return parameter values and optim output
  optim_output <- known_optim_output_stub(minustwolaploglik)

  # return list
  list(
    spcov_params_val = spcov_params_val, dispersion_params_val = dispersion_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, dispersion = dispersion_initial$is_known, randcov = randcov_initial$is_known)
  )
}
