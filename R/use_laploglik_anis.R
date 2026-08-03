# use_laploglik* family overview: see use_laploglik.R for the two axes that
# distinguish these sibling files (anisotropy, known-vs-estimated parameters).
# This file: anisotropic + estimated parameters, so an optim() search runs
# and, once finished, the two-candidate rotation heuristic described below is
# redone once more to determine which candidate optim() actually found.
#' Optimize the Laplace-approximated log-likelihood for GLM-type models under anisotropy
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param dispersion_initial A \code{dispersion_initial} object
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param spcov_profiled Whether the spatial covariance parameters are profiled out
#' @param randcov_initial A \code{randcov_initial} object (or \code{NULL} if there are no random effects)
#' @param randcov_profiled Whether the random effect variances are profiled out
#' @param optim_dotlist Additional arguments passed to \code{optim()}
#'
#' @return The same value as \code{use_laploglik()}, after also determining
#'   which of the rotation angle and its \code{pi}-complement optim()
#'   actually converged on (both are evaluated during the search and
#'   whichever fits better is kept -- an optimizer-robustness heuristic, not
#'   an identifiability correction; see the inline comment below)
#'
#' @noRd
use_laploglik_anis <- function(spcov_initial, dispersion_initial, data_object, estmethod, spcov_profiled,
                               randcov_initial = NULL, randcov_profiled = NULL, optim_dotlist) {
  # transforming to optim paramters (log odds or log scale)
  spcov_orig2optim_val <- spcov_orig2optim(
    spcov_initial = spcov_initial, spcov_profiled = spcov_profiled,
    data_object = data_object
  )

  # transforming to optim parameters
  dispersion_orig2optim_val <- dispersion_orig2optim(dispersion_initial)

  # transforming random effect parameters (if they are there else NULL)
  randcov_orig2optim_val <- randcov_orig2optim(
    randcov_initial = randcov_initial,
    randcov_profiled = randcov_profiled,
    spcov_initial = spcov_initial
  )

  # get optim par
  optim_par <- assemble_optim_par(
    spcov_orig2optim = spcov_orig2optim_val, randcov_orig2optim = randcov_orig2optim_val,
    dispersion_orig2optim = dispersion_orig2optim_val
  )

  # check optim dotlist
  optim_dotlist <- check_optim_method(optim_par, optim_dotlist)

  # performing optimization
  optim_output <- do.call("optim", c(
    list(
      par = optim_par,
      fn = laploglik_anis,
      spcov_orig2optim = spcov_orig2optim_val,
      dispersion_orig2optim = dispersion_orig2optim_val,
      data_object = data_object,
      estmethod = estmethod,
      spcov_profiled = spcov_profiled,
      randcov_orig2optim = randcov_orig2optim_val,
      randcov_profiled = randcov_profiled
    ),
    optim_dotlist
  ))

  # dispersion first then remove
  # the optimized par vector is packed as [spcov | randcov | dispersion] (see
  # assemble_optim_par()), so dispersion is peeled off the end first and
  # new_par (below) is what remains
  par <- optim_output$par
  dispersion_result <- peel_dispersion(dispersion_orig2optim_val, par, data_object$family)
  dispersion_params_val <- dispersion_result$dispersion_params_val
  par <- dispersion_result$par

  # reconcile a genuinely estimated ie with the numerical floor actually used to
  # build Sigma (see floor_estimated_ie()) -- done before the quadrant comparison
  # below so both candidate rotations are evaluated against the same ie
  unpacked <- unpack_optim2orig(spcov_orig2optim_val, randcov_orig2optim_val, par, spcov_profiled, randcov_profiled,
    data_object,
    spcov_initial = spcov_initial
  )
  spcov_params_val <- unpacked$spcov_params_val
  randcov_params_val <- unpacked$randcov_params_val

  resolved <- resolve_anis_rotation(spcov_params_val, randcov_params_val, data_object, estmethod,
    spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled,
    dispersion_params_val = dispersion_params_val
  )
  spcov_params_val <- resolved$spcov_params_val
  dist_matrix_list <- resolved$dist_matrix_list

  # return parameter values and optim output
  optim_output <- trim_optim_output(optim_output, optim_dotlist)

  # return list
  list(
    spcov_params_val = spcov_params_val, dispersion_params_val = dispersion_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, dispersion = dispersion_initial$is_known, randcov = randcov_initial$is_known)
  )
}
