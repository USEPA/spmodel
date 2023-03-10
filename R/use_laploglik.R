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
  optim_par <- get_optim_par_glm(spcov_orig2optim_val, dispersion_orig2optim_val, randcov_orig2optim_val)

  # check optim dotlist
  optim_dotlist <- check_optim_method(optim_par, optim_dotlist)

  # performing optimization
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
  par <- optim_output$par
  dispersion_orig_val <- dispersion_optim2orig(dispersion_orig2optim_val, par)
  dispersion_params_val <- dispersion_params(data_object$family, dispersion_orig_val$fill_orig_val)

  par <- dispersion_orig_val$new_par

  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim_val, par, spcov_profiled = spcov_profiled, data_object = data_object)

  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim_val), spcov_orig_val = spcov_orig_val)

  #
  # transforming to original scale
  randcov_orig_val <- randcov_optim2orig(randcov_orig2optim_val, spcov_orig2optim_val, par,
                                         randcov_profiled = randcov_profiled,
                                         spcov_optim2orig = spcov_params_val
  )

  # need to deal with list if randcov_profiled as sp variance changes
  # not used right now but could be
  # if (!is.null(randcov_profiled) && randcov_profiled) {
  #   spcov_params_val <- randcov_orig_val$spcov_optim2orig
  #   randcov_orig_val <- randcov_orig_val$fill_orig_val
  # }

  # making a random effects vector
  randcov_params_val <- randcov_params(randcov_orig_val)


  # # not used right now but could be
  # if (spcov_profiled && (is.null(randcov_profiled) ||
  #                        (!is.null(randcov_profiled) && randcov_profiled))) {
  #   # get the spcov_profiled variance
  #   sigma2 <- get_prof_sigma2(
  #     spcov_params_val, data_object, estmethod,
  #     dist_matrix_list, randcov_params_val
  #   )
  #
  #   # multiply by overall variance
  #   spcov_params_val[["de"]] <- sigma2 * spcov_params_val[["de"]]
  #   spcov_params_val[["ie"]] <- sigma2 * spcov_params_val[["ie"]]
  #
  #   if (!is.null(randcov_profiled)) {
  #     randcov_params_val <- sigma2 * randcov_params_val
  #   }
  #
  #   # add unconnected ar variance if needed
  #   if (inherits(spcov_params_val, c("car", "sar"))) {
  #     spcov_params_val[["extra"]] <- sigma2 * spcov_params_val[["extra"]]
  #   }
  # }

  # return parameter values and optim output
  optim_output <- list(
    method = optim_dotlist$method,
    control = optim_dotlist$control, value = optim_output$value,
    counts = optim_output$counts, convergence = optim_output$convergence,
    message = optim_output$message,
    hessian = if (optim_dotlist$hessian) optim_output$hessian else FALSE
  )
  # return list
  list(
    spcov_params_val = spcov_params_val, dispersion_params_val = dispersion_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, dispersion = dispersion_initial$is_known, randcov = randcov_initial$is_known)
  )
}
