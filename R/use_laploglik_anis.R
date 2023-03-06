use_laploglik_anis <- function(spcov_initial, dispersion_initial, data_object, estmethod, spcov_profiled,
                             randcov_initial = NULL, randcov_profiled = NULL, optim_dotlist) {

  # transforming to optim paramters (log odds or log scale)
  spcov_orig2optim_val <- spcov_orig2optim(spcov_initial = spcov_initial, spcov_profiled = spcov_profiled)

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

}
