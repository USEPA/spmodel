laploglik_anis <- function(par, spcov_orig2optim, dispersion_orig2optim, data_object, estmethod,
                         spcov_profiled, randcov_orig2optim = NULL,
                         randcov_profiled = NULL) {

  # dispersion first then remove
  dispersion_orig_val <- dispersion_optim2orig(dispersion_orig2optim, par)
  dispersion_params_val <- dispersion_params(data_object$family, dispersion_orig_val$fill_orig_val)

  par <- dispersion_orig_val$new_par

  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim, par, spcov_profiled = spcov_profiled)

  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim), spcov_orig_val = spcov_orig_val)

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

  # quadrant 1
  ## make distance matrix
  new_coords_list_q1 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
                               rotate = spcov_params_val[["rotate"]], scale = spcov_params_val[["scale"]]
  )

  dist_matrix_list_q1 <- lapply(new_coords_list_q1, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

  # compute relevant products
  lapll_prods_q1 <- laploglik_products(
    spcov_params_val, dispersion_params_val, data_object, estmethod,
    dist_matrix_list_q1, randcov_params_val
  )

  # find -2loglik
  minustwolaploglik_q1 <- get_minustwolaploglik(lapll_prods_q1, estmethod, data_object$n, data_object$p, spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled)


  new_coords_list_q2 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
                               rotate = abs(pi - spcov_params_val[["rotate"]]), scale = spcov_params_val[["scale"]]
  )
  dist_matrix_list_q2 <- lapply(new_coords_list_q2, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

  # compute relevant products
  lapll_prods_q2 <- laploglik_products(
    spcov_params_val, dispersion_params_val, data_object, estmethod,
    dist_matrix_list_q2, randcov_params_val
  )

  # find -2loglik
  minustwolaploglik_q2 <- get_minustwolaploglik(lapll_prods_q2, estmethod, data_object$n,
                                          data_object$p,
                                          spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled
  )

  minustwolaploglik <- min(c(minustwolaploglik_q1, minustwolaploglik_q2))
}
