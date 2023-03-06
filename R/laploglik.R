laploglik <- function(par, spcov_orig2optim, dispersion_orig2optim, data_object, estmethod, dist_matrix_list,
                    spcov_profiled, randcov_orig2optim = NULL,
                    randcov_profiled = NULL) {


  # dispersion first then remove
  dispersion_orig_val <- dispersion_optim2orig(dispersion_orig2optim, par)
  dispersion_params_val <- dispersion_params(data_object$family, dispersion_orig_val)

  par <- par[-length(par)] # remove dispersion par

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
  lapll_prods <- laploglik_products(
    spcov_params_val, dispersion_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )

  minustwolaploglik <- get_minustwolaploglik(lapll_prods, estmethod, data_object$n,
                                          data_object$p, spcov_profiled = FALSE, randcov_profiled = FALSE)
}
