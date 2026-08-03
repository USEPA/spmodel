#' Gaussian log-likelihood estimation for geostatistical models
#'
#' @param data_object The data object
#' @param formula A formula
#' @param spcov_initial The spatial initial object
#' @param estmethod The estimation method
#' @param optim_dotlist The optim dotlist
#'
#' @return The Gaussian log-likelihood estimates
#'
#' @noRd
cov_estimate_laploglik_spglm <- function(data_object, formula, spcov_initial,
                                         dispersion_initial, estmethod, optim_dotlist) {
  # For non-Gaussian responses the exact marginal likelihood generally has no closed
  # form (it requires integrating out spatial random effects), so a Laplace
  # approximation to the log-likelihood is used instead; this mirrors
  # cov_estimate_gloglik_splm() but branches only on "known" vs. "not known" since
  # there is no cheap iid special case for the Laplace approximation
  # make NA spcov_initial
  spcov_initial_NA_val <- spcov_initial_NA_glm(data_object$family, spcov_initial, anisotropy = data_object$anisotropy)

  # make NA dispersion initial
  dispersion_initial_NA_val <- dispersion_initial_NA(dispersion_initial, data_object)

  # store distance matrix (if applicable)
  if (data_object$anisotropy) {
    dist_matrix_list <- NULL
  } else {
    dist_matrix_list <- lapply(data_object$obdata_list, function(x) spdist(x, data_object$xcoord, data_object$ycoord))
  }

  if (is.null(data_object$randcov_initial)) {
    cov_initial_val <- cov_initial_search_glm(
      spcov_initial_NA = spcov_initial_NA_val,
      dispersion_initial_NA = dispersion_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list
    )

    # spatial and random effect initial values
    spcov_initial_val <- cov_initial_val$spcov_initial_val
    dispersion_initial_val <- cov_initial_val$dispersion_initial_val
    randcov_initial_val <- NULL
  } else {
    randcov_names <- data_object$randcov_names
    randcov_initial_NA_val <- randcov_initial_NA(data_object$randcov_initial, randcov_names)


    cov_initial_val <- cov_initial_search_glm(
      spcov_initial_NA = spcov_initial_NA_val,
      dispersion_initial_NA = dispersion_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list,
      randcov_initial_NA = randcov_initial_NA_val
    )

    # spatial and random effect initial values
    spcov_initial_val <- cov_initial_val$spcov_initial_val
    dispersion_initial_val <- cov_initial_val$dispersion_initial_val
    randcov_initial_val <- cov_initial_val$randcov_initial_val
  }

  # choose among known/non-profiled likelihood evaluators -- see
  # run_laploglik_dispatch_spglm() in cov_estimate_dispatch_helpers.R
  cov_estimate_val <- run_laploglik_dispatch_spglm(
    spcov_initial_val, dispersion_initial_val, randcov_initial_val, data_object, estmethod, dist_matrix_list, optim_dotlist
  )
}

#' Gaussian log-likelihood estimation for areal (autoregressive) GLM models
#'
#' @param data_object The data object
#' @param formula A formula
#' @param spcov_initial The spatial initial object
#' @param dispersion_initial The dispersion initial object
#' @param estmethod The estimation method
#' @param optim_dotlist The optim dotlist
#'
#' @return The Gaussian log-likelihood estimates
#'
#' @noRd
cov_estimate_laploglik_spgautor <- function(data_object, formula, spcov_initial,
                                            dispersion_initial, estmethod,
                                            optim_dotlist) {
  # areal (CAR/SAR) counterpart to cov_estimate_laploglik_spglm(); uses the
  # neighbor weight matrix W in place of a distance matrix
  # make NA spcov_initial
  spcov_initial_NA_val <- spcov_initial_NA_glm(data_object$family, spcov_initial, is_W_connected = data_object$is_W_connected)

  # make NA dispersion initial
  dispersion_initial_NA_val <- dispersion_initial_NA(dispersion_initial, data_object)

  # make dist_matrix_list
  # NOTE THIS IS NOT ACTUALLY A LIST WITH SPAUTO() BUT NAME
  # KEPT FOR CONSISTENCY WITH SPLM()
  dist_matrix_list <- data_object$W

  if (is.null(data_object$randcov_initial)) {
    # find initial values
    cov_initial_val <- cov_initial_search_glm(
      spcov_initial_NA = spcov_initial_NA_val,
      dispersion_initial_NA = dispersion_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list
    )

    # spatial and random effect initial values
    spcov_initial_val <- cov_initial_val$spcov_initial_val
    dispersion_initial_val <- cov_initial_val$dispersion_initial_val
    randcov_initial_val <- cov_initial_val$randcov_initial_val
  } else {
    # assign random effects
    randcov_names <- data_object$randcov_names
    randcov_initial_NA_val <- randcov_initial_NA(data_object$randcov_initial, randcov_names)

    # find initial values
    cov_initial_val <- cov_initial_search_glm(
      spcov_initial_NA = spcov_initial_NA_val,
      dispersion_initial_NA = dispersion_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list,
      randcov_initial_NA = randcov_initial_NA_val
    )

    # spatial and random effect initial values
    spcov_initial_val <- cov_initial_val$spcov_initial_val
    dispersion_initial_val <- cov_initial_val$dispersion_initial_val
    randcov_initial_val <- cov_initial_val$randcov_initial_val
  }

  # choose among known/non-profiled likelihood evaluators -- see
  # run_laploglik_dispatch_spgautor() in cov_estimate_dispatch_helpers.R
  cov_estimate_val <- run_laploglik_dispatch_spgautor(
    spcov_initial_val, dispersion_initial_val, randcov_initial_val, data_object, estmethod, dist_matrix_list, optim_dotlist
  )
}
