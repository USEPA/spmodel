use_laploglik_known <- function(spcov_initial, dispersion_initial, data_object, estmethod, dist_matrix_list, randcov_initial) {

  spcov_params_val <- get_spcov_params(class(spcov_initial), spcov_initial$initial)
  dispersion_params_val <- dispersion_params(data_object$family, unname(dispersion_initial$initial))
  # unname otherwise name dispersion.dispersion
  randcov_params_val <- randcov_params(randcov_initial$initial)

  lapll_prods <- laploglik_products(
    spcov_params_val, dispersion_params_val, data_object, estmethod,
    dist_matrix_list, randcov_params_val
  )
  ## compute -2ll
  minustwolaploglik <- get_minustwolaploglik(lapll_prods, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE)
  # return parameter values and optim output
  optim_output <- list(
    method = NA, control = NA, value = minustwolaploglik,
    counts = NA, convergence = NA,
    message = NA, hessian = NA
  )



  # return list
  list(
    spcov_params_val = spcov_params_val, dispersion_params_val = dispersion_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, dispersion = dispersion_initial$is_known, randcov = randcov_initial$is_known)
  )

}
