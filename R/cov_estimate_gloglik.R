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

  # make NA spcov_initial
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy = data_object$anisotropy)

  # store distance matrix (if applicable)
  if (data_object$anisotropy) {
    dist_matrix_list <- NULL
  } else {
    if (inherits(spcov_initial, "none") && is.null(data_object$randcov_initial)) {
      dist_matrix_list <- NULL
    } else {
      dist_matrix_list <- lapply(data_object$obdata_list, function(x) spdist(x, data_object$xcoord, data_object$ycoord))
    }
  }

  if (is.null(data_object$randcov_initial)) {
    cov_initial_val <- cov_initial_search(
      spcov_initial_NA = spcov_initial_NA_val,
      estmethod = estmethod,
      data_object = data_object,
      dist_matrix_list = dist_matrix_list
    )

    spcov_initial_val <- cov_initial_val$spcov_initial_val

    # some initial logical statements to figure out what likelihood to use
    de_known <- spcov_initial_val$is_known[["de"]]
    de_known_zero <- de_known && (spcov_initial_val$initial[["de"]] == 0)
    ie_known <- spcov_initial_val$is_known[["ie"]]

    # all parameters known
    if (all(spcov_initial_val$is_known)) {
      if (data_object$anisotropy) {
        cov_estimate_val <- use_gloglik_known_anis(spcov_initial_val, data_object, estmethod, randcov_initial = NULL)
      } else {
        cov_estimate_val <- use_gloglik_known(spcov_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial = NULL)
      }
    } else if (de_known_zero) { # de = 0 for iid models
      cov_estimate_val <- use_gloglik_iid(spcov_initial_val, estmethod, data_object, dist_matrix_list)
    } else if (de_known || ie_known) { # use non-profiled likelihood
      if (data_object$anisotropy) {
        cov_estimate_val <- use_gloglik_anis(spcov_initial_val, data_object, estmethod,
          spcov_profiled = FALSE, optim_dotlist = optim_dotlist
        )
      } else {
        cov_estimate_val <- use_gloglik(spcov_initial_val, data_object, estmethod,
          dist_matrix_list,
          spcov_profiled = FALSE, optim_dotlist = optim_dotlist
        )
      }
    } else {
      if (data_object$anisotropy) {
        cov_estimate_val <- use_gloglik_anis(spcov_initial_val, data_object, estmethod,
          spcov_profiled = TRUE, optim_dotlist = optim_dotlist
        )
      } else {
        cov_estimate_val <- use_gloglik(spcov_initial_val, data_object, estmethod,
          dist_matrix_list,
          spcov_profiled = TRUE, optim_dotlist = optim_dotlist
        )
      }
    }
  } else {
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

    # some initial logical statements to figure out what likelihood to use
    de_known <- spcov_initial_val$is_known[["de"]]
    de_known_zero <- de_known && (spcov_initial_val$initial[["de"]] == 0)
    ie_known <- spcov_initial_val$is_known[["ie"]]

    if (all(spcov_initial_val$is_known, randcov_initial_val$is_known)) {
      if (data_object$anisotropy) {
        cov_estimate_val <- use_gloglik_known_anis(spcov_initial_val, data_object, estmethod, randcov_initial_val)
      } else {
        cov_estimate_val <- use_gloglik_known(spcov_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial_val)
      }
    } else if (any(de_known && !de_known_zero, spcov_initial_val$is_known[["ie"]], randcov_initial_val$is_known)) { # at least one variance parameter known
      if (data_object$anisotropy) {
        cov_estimate_val <- use_gloglik_anis(spcov_initial_val, data_object, estmethod,
          spcov_profiled = FALSE,
          randcov_initial = randcov_initial_val, randcov_profiled = FALSE,
          optim_dotlist = optim_dotlist
        )
      } else {
        cov_estimate_val <- use_gloglik(spcov_initial_val, data_object, estmethod,
          dist_matrix_list = dist_matrix_list, spcov_profiled = FALSE,
          randcov_initial = randcov_initial_val, randcov_profiled = FALSE,
          optim_dotlist = optim_dotlist
        )
      }
    } else { # all unknown so profile
      if (data_object$anisotropy) {
        cov_estimate_val <- use_gloglik_anis(spcov_initial_val, data_object, estmethod,
          spcov_profiled = TRUE,
          randcov_initial = randcov_initial_val, randcov_profiled = TRUE,
          optim_dotlist = optim_dotlist
        )
      } else {
        cov_estimate_val <- use_gloglik(spcov_initial_val, data_object, estmethod,
          dist_matrix_list = dist_matrix_list, spcov_profiled = TRUE,
          randcov_initial = randcov_initial_val, randcov_profiled = TRUE,
          optim_dotlist = optim_dotlist
        )
      }
    }
  }
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

    # some initial logical statements to figure out what likelihood to use
    de_known <- spcov_initial_val$is_known[["de"]]
    de_known_zero <- de_known && (spcov_initial_val$initial[["de"]] == 0)
    ie_known <- spcov_initial_val$is_known[["ie"]]
    extra_known <- spcov_initial_val$is_known[["extra"]]
    extra_known_zero <- extra_known && (spcov_initial_val$initial[["extra"]] == 0)

    # use likelihood
    if (all(spcov_initial_val$is_known)) {
      cov_estimate_val <- use_gloglik_known(spcov_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial = NULL)
    } else if (de_known_zero && extra_known_zero) {
      cov_estimate_val <- use_gloglik_iid(spcov_initial_val, estmethod, data_object, dist_matrix_list)
    } else if (extra_known_zero && !de_known && !ie_known) {
      cov_estimate_val <- use_gloglik(spcov_initial_val, data_object, estmethod,
        dist_matrix_list,
        spcov_profiled = TRUE,
        optim_dotlist = optim_dotlist
      )
    } else {
      cov_estimate_val <- use_gloglik(spcov_initial_val, data_object, estmethod,
        dist_matrix_list,
        spcov_profiled = FALSE,
        optim_dotlist = optim_dotlist
      )
    }
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

    # some initial logical statements to figure out what likelihood to use
    de_known <- spcov_initial_val$is_known[["de"]]
    de_known_zero <- de_known && (spcov_initial_val$initial[["de"]] == 0)
    ie_known <- spcov_initial_val$is_known[["ie"]]
    extra_known <- spcov_initial_val$is_known[["extra"]]
    extra_known_zero <- extra_known && (spcov_initial_val$initial[["extra"]] == 0)

    # use likelihood
    if (all(spcov_initial_val$is_known, randcov_initial_val$is_known)) {
      cov_estimate_val <- use_gloglik_known(spcov_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial_val)
    } else if (extra_known_zero && !any(de_known, ie_known, randcov_initial_val$is_known)) {
      cov_estimate_val <- use_gloglik(spcov_initial_val, data_object, estmethod,
        dist_matrix_list = dist_matrix_list, spcov_profiled = TRUE,
        randcov_initial = randcov_initial_val, randcov_profiled = TRUE,
        optim_dotlist = optim_dotlist
      )
    } else {
      cov_estimate_val <- use_gloglik(spcov_initial_val, data_object, estmethod,
        dist_matrix_list = dist_matrix_list, spcov_profiled = FALSE,
        randcov_initial = randcov_initial_val, randcov_profiled = FALSE,
        optim_dotlist = optim_dotlist
      )
    }
  }
}
