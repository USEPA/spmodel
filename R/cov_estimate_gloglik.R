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
cov_estimate_gloglik_splm <- function(data_object, formula, spcov_initial, estmethod,
                                      optim_dotlist) {
  # This function chooses among several specialized Gaussian log-likelihood
  # evaluators (known, iid, profiled, non-profiled, anisotropic) rather than always
  # calling one general routine, because each specialized path skips computation
  # that isn't needed for that case (e.g. no optimization at all if every parameter
  # is fixed, or a simpler iid likelihood when the spatial dependence is zero).

  # mark which spatial covariance parameters the user fixed (vs. left to estimate)
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy = data_object$anisotropy)

  # store distance matrix (if applicable); not needed at all for "none"/"ie"
  # covariance types with no random effects, since those models have no spatial
  # dependence structure that depends on distance
  if (data_object$anisotropy) {
    dist_matrix_list <- NULL
  } else {
    if (inherits(spcov_initial, c("none", "ie")) && is.null(data_object$randcov_initial)) {
      dist_matrix_list <- NULL
    } else {
      dist_matrix_list <- lapply(data_object$obdata_list, function(x) spdist(x, data_object$xcoord, data_object$ycoord))
    }
  }

  if (is.null(data_object$randcov_initial)) {
    # grid search for good optimizer starting values
    cov_initial_val <- cov_initial_search(
      spcov_initial_NA = spcov_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list
    )

    spcov_initial_val <- cov_initial_val$spcov_initial_val
    randcov_initial_val <- NULL
  } else {
    # branch for models that include one or more random effects in addition to the
    # spatial covariance structure
    randcov_names <- data_object$randcov_names
    randcov_initial_NA_val <- randcov_initial_NA(data_object$randcov_initial, randcov_names)

    cov_initial_val <- cov_initial_search(
      spcov_initial_NA = spcov_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list,
      randcov_initial_NA = randcov_initial_NA_val
    )

    spcov_initial_val <- cov_initial_val$spcov_initial_val
    randcov_initial_val <- cov_initial_val$randcov_initial_val
  }

  # choose among known/iid/non-profiled/profiled likelihood evaluators -- see
  # run_gloglik_dispatch_splm() in cov_estimate_dispatch_helpers.R
  cov_estimate_val <- run_gloglik_dispatch_splm(
    spcov_initial_val, randcov_initial_val, data_object, estmethod, dist_matrix_list, optim_dotlist
  )
}


#' Gaussian log-likelihood estimation for autoregressive models
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
cov_estimate_gloglik_spautor <- function(data_object, formula, spcov_initial, estmethod,
                                         optim_dotlist) {
  # For autoregressive (CAR/SAR) models the "extra" parameter is the autocorrelation
  # parameter (rho); extra = 0 means no spatial dependence, in which case (combined
  # with de = 0) the model reduces to the plain iid case, same as for geostatistical
  # models
  # make NA spcov_initial
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, is_W_connected = data_object$is_W_connected)

  # make dist_matrix_list
  # NOTE THIS IS NOT ACTUALLY A LIST WITH SPAUTO() BUT NAME
  # KEPT FOR CONSISTENCY WITH SPLM()
  dist_matrix_list <- data_object$W

  if (is.null(data_object$randcov_initial)) {
    # find initial values
    cov_initial_val <- cov_initial_search(
      spcov_initial_NA = spcov_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list
    )

    # initial spatial covariance value
    spcov_initial_val <- cov_initial_val$spcov_initial_val
    randcov_initial_val <- NULL
  } else {
    # assign random effects
    randcov_names <- data_object$randcov_names
    randcov_initial_NA_val <- randcov_initial_NA(data_object$randcov_initial, randcov_names)

    # find initial values
    cov_initial_val <- cov_initial_search(
      spcov_initial_NA = spcov_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list,
      randcov_initial_NA = randcov_initial_NA_val
    )

    # spatial and random effect initial values
    spcov_initial_val <- cov_initial_val$spcov_initial_val
    randcov_initial_val <- cov_initial_val$randcov_initial_val
  }

  # choose among known/iid/profiled/non-profiled likelihood evaluators -- see
  # run_gloglik_dispatch_spautor() in cov_estimate_dispatch_helpers.R
  cov_estimate_val <- run_gloglik_dispatch_spautor(
    spcov_initial_val, randcov_initial_val, data_object, estmethod, dist_matrix_list, optim_dotlist
  )
}
