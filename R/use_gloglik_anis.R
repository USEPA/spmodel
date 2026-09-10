#' Use Gaussian log-likelihood estimation with anisotropy
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param estmethod The estimation method (\code{"reml"} or \code{"ml"})
#' @param X Model matrix
#' @param y Response vector
#' @param n Sample size
#' @param p Number of fixed effects
#' @param xcoord_val x-coordinate value
#' @param ycoord_val y-coordinate value
#' @param spcov_profiled Is the spatial covariance profiled?
#' @param optim_dotlist optim dotlist
#' @param randcov_initial A \code{randcov_initial} object
#' @param randcov_Zs Random effects design matrices
#' @param observed_index The index of observed values
#' @param partition_matrix The partition matrix
#'
#' @return Estimated covariance parameters
#'
#' @noRd
use_gloglik_anis <- function(spcov_initial, data_object, estmethod, spcov_profiled,
                             randcov_initial = NULL, randcov_profiled = NULL, optim_dotlist) {
  # transforming to optim paramters (log odds or log scale)
  spcov_orig2optim_val <- spcov_orig2optim(
    spcov_initial = spcov_initial, spcov_profiled = spcov_profiled,
    data_object = data_object
  )


  # transforming random effect parameters (if they are there else NULL)
  randcov_orig2optim_val <- randcov_orig2optim(
    randcov_initial = randcov_initial,
    randcov_profiled = randcov_profiled,
    spcov_initial = spcov_initial
  )


  # get optim par
  optim_par <- assemble_optim_par(spcov_orig2optim = spcov_orig2optim_val, randcov_orig2optim = randcov_orig2optim_val)

  # check optim dotlist
  optim_dotlist <- check_optim_method(optim_par, optim_dotlist)

  # performing optimization
  optim_output <- do.call("optim", c(
    list(
      par = optim_par,
      fn = gloglik_anis,
      spcov_orig2optim = spcov_orig2optim_val,
      data_object = data_object,
      estmethod = estmethod,
      spcov_profiled = spcov_profiled,
      randcov_orig2optim = randcov_orig2optim_val,
      randcov_profiled = randcov_profiled
    ),
    optim_dotlist
  ))

  # transforming to original scale. NOTE: unlike the laploglik-family call
  # sites, spcov_initial is deliberately NOT passed here -- floor_estimated_ie()
  # must run after the profiled-variance rescale below (diagtol is an
  # absolute-scale threshold, and ie is not on its final absolute scale until
  # after that rescale), so it is called separately further down
  unpacked <- unpack_optim2orig(spcov_orig2optim_val, randcov_orig2optim_val, optim_output$par, spcov_profiled, randcov_profiled, data_object)
  spcov_params_val <- unpacked$spcov_params_val
  randcov_params_val <- unpacked$randcov_params_val

  resolved <- resolve_anis_rotation(spcov_params_val, randcov_params_val, data_object, estmethod,
    spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled
  )
  spcov_params_val <- resolved$spcov_params_val
  dist_matrix_list <- resolved$dist_matrix_list

  rescaled <- rescale_profiled_variance(spcov_params_val, randcov_params_val, data_object, estmethod,
    dist_matrix_list, spcov_profiled, randcov_profiled
  )
  spcov_params_val <- rescaled$spcov_params_val
  randcov_params_val <- rescaled$randcov_params_val

  # reconcile a genuinely estimated ie with the numerical floor
  # spcov_matrix.*() applies internally when building Sigma -- mirrors the
  # equivalent GLM-side reconciliation, see R/floor_estimated_ie.R. Must run
  # here, after the profiled-variance rescale above, not inside
  # unpack_optim2orig() like the laploglik call sites do: diagtol is an
  # absolute-variance-scale threshold, but ie is still a correlation-scale
  # proportion (summing to 1 with de) until the rescale multiplies it back
  # into absolute units.
  spcov_params_val <- floor_estimated_ie(spcov_params_val, spcov_initial$is_known, data_object$diagtol)

  # return parameter values and optim output
  optim_output <- trim_optim_output(optim_output, optim_dotlist)


  # return list
  list(
    spcov_params_val = spcov_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, randcov = randcov_initial$is_known)
  )
}
