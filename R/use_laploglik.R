# use_laploglik* family overview: these functions optimize (or, for "_known"
# variants, simply evaluate) the Laplace-approximated log-likelihood used to
# fit GLM-type spatial models (spglm/spgautor). Like the use_gloglik* family,
# variants differ along anisotropy (plain vs "_anis") and known-vs-estimated
# parameters (plain vs "_known"). This file: isotropic + estimated parameters,
# so a full optim() search runs over spatial covariance, dispersion, and
# random effect parameters jointly. Compare with use_laploglik_anis.R
# (anisotropic, estimated), use_laploglik_known.R (isotropic, known), and
# use_laploglik_known_anis.R (anisotropic, known).
#' Optimize the Laplace-approximated log-likelihood for GLM-type models
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param dispersion_initial A \code{dispersion_initial} object
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param dist_matrix_list A list of distance matrices
#' @param spcov_profiled Whether the spatial covariance parameters are profiled out
#' @param randcov_initial A \code{randcov_initial} object (or \code{NULL} if there are no random effects)
#' @param randcov_profiled Whether the random effect variances are profiled out
#' @param optim_dotlist Additional arguments passed to \code{optim()}
#'
#' @return A list with the estimated spatial covariance, dispersion, and
#'   random effect parameters; the \code{optim()} output; the distance matrix
#'   list; and which parameters were fixed (\code{is_known})
#'
#' @noRd
use_laploglik <- function(spcov_initial, dispersion_initial, data_object, estmethod, dist_matrix_list, spcov_profiled,
                          randcov_initial = NULL, randcov_profiled = NULL, optim_dotlist) {
  # transforming to optim paramters (log odds or log scale)
  spcov_orig2optim_val <- spcov_orig2optim(spcov_initial = spcov_initial, spcov_profiled = spcov_profiled, data_object = data_object)

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
  # optim_par packs the spcov, dispersion, and randcov parameters (all on
  # their unconstrained optim scales) into one vector that laploglik()
  # unpacks and jointly optimizes
  optim_output <- do.call("optim", c(
    list(
      par = optim_par,
      fn = laploglik,
      spcov_orig2optim = spcov_orig2optim_val,
      dispersion_orig2optim = dispersion_orig2optim_val,
      data_object = data_object,
      estmethod = estmethod,
      dist_matrix_list = dist_matrix_list,
      spcov_profiled = spcov_profiled,
      randcov_orig2optim = randcov_orig2optim_val,
      randcov_profiled = randcov_profiled
    ),
    optim_dotlist
  ))

  # dispersion first then remove
  # the optimized par vector is packed as [spcov | randcov | dispersion] (see
  # assemble_optim_par()), so dispersion is peeled off the end (converted back
  # to its original scale) first, and new_par is what remains (spcov +
  # randcov) for the steps below
  par <- optim_output$par
  dispersion_result <- peel_dispersion(dispersion_orig2optim_val, par, data_object$family)
  dispersion_params_val <- dispersion_result$dispersion_params_val
  par <- dispersion_result$par

  unpacked <- unpack_optim2orig(spcov_orig2optim_val, randcov_orig2optim_val, par, spcov_profiled, randcov_profiled,
    data_object,
    spcov_initial = spcov_initial
  )
  spcov_params_val <- unpacked$spcov_params_val
  randcov_params_val <- unpacked$randcov_params_val

  # return parameter values and optim output
  optim_output <- trim_optim_output(optim_output, optim_dotlist)
  # return list
  list(
    spcov_params_val = spcov_params_val, dispersion_params_val = dispersion_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, dispersion = dispersion_initial$is_known, randcov = randcov_initial$is_known)
  )
}
