#' Find (minus twice negative) Gaussian log-likelihood while optimizing (with anisotropy)
#'
#' @param par Parameters to optimize over
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param estmethod The estimation method
#' @param X Model matrix
#' @param y Response vector
#' @param n Sample size
#' @param p Number of fixed effects
#' @param xcoord_val A vector of x-coordinates
#' @param ycoord_val A vector of y-coordinates
#' @param spcov_profiled Whether the overall spatial variance is profiled
#' @param randcov_orig2optim A \code{randcov_orig2optim} object
#' @param randcov_Zs Random effect design matrices
#' @param observed_index Index of observed values
#' @param partition_matrix Partition matrix
#'
#' @return (Minus twice negative) Gaussian log-likelihood (with anisotropy)
#'
#' @noRd
gloglik_anis <- function(par, spcov_orig2optim, data_object, estmethod,
                         spcov_profiled, randcov_orig2optim = NULL,
                         randcov_profiled = NULL) {

  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim, par, spcov_profiled = spcov_profiled,
                                     data_object = data_object)

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

  minustwologlik <- min(c(minustwologlik_q1, minustwologlik_q2))
}
