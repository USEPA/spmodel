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
  spcov_orig2optim_val <- spcov_orig2optim(spcov_initial = spcov_initial, spcov_profiled = spcov_profiled)


  # browser()
  # transforming random effect parameters (if they are there else NULL)
  randcov_orig2optim_val <- randcov_orig2optim(
    randcov_initial = randcov_initial,
    randcov_profiled = randcov_profiled,
    spcov_initial = spcov_initial
  )


  # get optim par
  optim_par <- get_optim_par(spcov_orig2optim_val, randcov_orig2optim_val)

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

  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim_val, optim_output$par, spcov_profiled = spcov_profiled)

  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim_val), spcov_orig_val = spcov_orig_val)

  # transforming to original scale
  randcov_orig_val <- randcov_optim2orig(randcov_orig2optim_val,
    spcov_orig2optim_val,
    optim_output$par,
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

  # finding appropriate rotation
  # quadrant 1

  ## make distance matrix
  new_coords_list_q1 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
    rotate = spcov_params_val[["rotate"]], scale = spcov_params_val[["scale"]]
  )

  dist_matrix_list_q1 <- lapply(new_coords_list_q1, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

  # compute relevant products
  gll_prods_q1 <- gloglik_products(
    spcov_params_val, data_object, estmethod,
    dist_matrix_list_q1, randcov_params_val
  )

  # find -2loglik
  minustwologlik_q1 <- get_minustwologlik(gll_prods_q1, estmethod, data_object$n, data_object$p, spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled)


  new_coords_list_q2 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
    rotate = abs(pi - spcov_params_val[["rotate"]]), scale = spcov_params_val[["scale"]]
  )
  dist_matrix_list_q2 <- lapply(new_coords_list_q2, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

  # compute relevant products
  gll_prods_q2 <- gloglik_products(
    spcov_params_val, data_object, estmethod,
    dist_matrix_list_q2, randcov_params_val
  )

  # find -2loglik
  minustwologlik_q2 <- get_minustwologlik(gll_prods_q2, estmethod, data_object$n,
    data_object$p,
    spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled
  )

  ## find appropriate value
  rotate_min <- which.min(c(minustwologlik_q1, minustwologlik_q2))

  if (rotate_min == 1) {
    dist_matrix_list <- dist_matrix_list_q1
  } else if (rotate_min == 2) {
    spcov_params_val[["rotate"]] <- abs(pi - spcov_params_val[["rotate"]])
    dist_matrix_list <- dist_matrix_list_q2
  }

  # browser()
  if (spcov_profiled && (is.null(randcov_profiled) ||
    (!is.null(randcov_profiled) && randcov_profiled))) {
    # get the spcov_profiled variance
    sigma2 <- get_prof_sigma2(
      spcov_params_val, data_object, estmethod,
      dist_matrix_list, randcov_params_val
    )

    # multiply by overall variance
    spcov_params_val[["de"]] <- sigma2 * spcov_params_val[["de"]]
    spcov_params_val[["ie"]] <- sigma2 * spcov_params_val[["ie"]]

    if (!is.null(randcov_profiled)) {
      randcov_params_val <- sigma2 * randcov_params_val
    }

    # add unconnected ar variance if needed
    if (inherits(spcov_params_val, c("car", "sar"))) {
      spcov_params_val[["extra"]] <- sigma2 * spcov_params_val[["extra"]]
    }
  }

  # return parameter values and optim output
  optim_output <- list(
    method = optim_dotlist$method, control = optim_dotlist$control,
    value = optim_output$value,
    counts = optim_output$counts, convergence = optim_output$convergence,
    message = optim_output$message, hessian = if (optim_dotlist$hessian) optim_output$hessian else FALSE
  )


  # return list
  list(
    spcov_params_val = spcov_params_val, randcov_params_val = randcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known, randcov = randcov_initial$is_known)
  )
}
