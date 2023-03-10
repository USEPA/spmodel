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


  # make NA spcov_initial
  spcov_initial_NA_val <- spcov_initial_NA_glm(data_object$family, spcov_initial, anisotropy = data_object$anisotropy)
  # dispersion ie confounded in comment below
  # spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy = data_object$anisotropy)

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

    # some initial logical statements to figure out what likelihood to use
    de_known <- spcov_initial_val$is_known[["de"]]
    ie_known <- spcov_initial_val$is_known[["ie"]]
    dispersion_known <- dispersion_initial_val$is_known[["dispersion"]]

    if (all(de_known, ie_known, dispersion_known)) {
      if (data_object$anisotropy) {
        cov_estimate_val <- use_laploglik_known_anis(spcov_initial_val, dispersion_initial_val, data_object, estmethod, randcov_initial = NULL)
      } else {
        cov_estimate_val <- use_laploglik_known(spcov_initial_val, dispersion_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial = NULL)
      }
    } else {
      if (data_object$anisotropy) {
        cov_estimate_val <- use_laploglik_anis(spcov_initial_val, dispersion_initial_val, data_object, estmethod,
                                             spcov_profiled = FALSE, optim_dotlist = optim_dotlist
        )
      } else {
        cov_estimate_val <- use_laploglik(spcov_initial_val, dispersion_initial_val, data_object, estmethod,
                                        dist_matrix_list,
                                        spcov_profiled = FALSE, optim_dotlist = optim_dotlist
        )
      }
    }
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

    # some initial logical statements to figure out what likelihood to use
    de_known <- spcov_initial_val$is_known[["de"]]
    ie_known <- spcov_initial_val$is_known[["ie"]]
    dispersion_known <- dispersion_initial_val$is_known[["dispersion"]]

    if (all(de_known, ie_known, dispersion_known, randcov_initial_val$is_known)) {
      if (data_object$anisotropy) {
        cov_estimate_val <- use_laploglik_known_anis(spcov_initial_val, dispersion_initial_val, data_object, estmethod, randcov_initial_val)
      } else {
        cov_estimate_val <- use_laploglik_known(spcov_initial_val, dispersion_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial_val)
      }
    } else {
      if (data_object$anisotropy) {
        cov_estimate_val <- use_laploglik_anis(spcov_initial_val, dispersion_initial_val, data_object, estmethod,
                                             spcov_profiled = FALSE,
                                             randcov_initial = randcov_initial_val, randcov_profiled = FALSE,
                                             optim_dotlist = optim_dotlist
        )
      } else {
        cov_estimate_val <- use_laploglik(spcov_initial_val, dispersion_initial_val, data_object, estmethod,
                                        dist_matrix_list = dist_matrix_list, spcov_profiled = FALSE,
                                        randcov_initial = randcov_initial_val, randcov_profiled = FALSE,
                                        optim_dotlist = optim_dotlist
        )
      }
    }
  }
}

cov_estimate_laploglik_spgautor <- function(data_object, formula, spcov_initial,
                                         dispersion_initial, estmethod,
                                         optim_dotlist) {
  # make NA spcov_initial
  spcov_initial_NA_val <- spcov_initial_NA_glm(data_object$family, spcov_initial, is_W_connected = data_object$is_W_connected)
  # dispersion ie confounded in comment below
  # spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, is_W_connected = data_object$is_W_connected)

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

    # some initial logical statements to figure out what likelihood to use
    de_known <- spcov_initial_val$is_known[["de"]]
    ie_known <- spcov_initial_val$is_known[["ie"]]
    extra_known <- spcov_initial_val$is_known[["extra"]]
    dispersion_known <- dispersion_initial_val$is_known[["dispersion"]]

    if (all(de_known, ie_known, extra_known, dispersion_known)) {
      cov_estimate_val <- use_laploglik_known(spcov_initial_val, dispersion_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial = NULL)
    } else {
      cov_estimate_val <- use_laploglik(spcov_initial_val, dispersion_initial_val, data_object, estmethod,
                                        dist_matrix_list,
                                        spcov_profiled = FALSE, optim_dotlist = optim_dotlist
      )
    }
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

    # some initial logical statements to figure out what likelihood to use
    de_known <- spcov_initial_val$is_known[["de"]]
    ie_known <- spcov_initial_val$is_known[["ie"]]
    extra_known <- spcov_initial_val$is_known[["extra"]]
    dispersion_known <- dispersion_initial_val$is_known[["dispersion"]]

    if (all(de_known, ie_known, extra_known, dispersion_known, randcov_initial_val$is_known)) {
      cov_estimate_val <- use_laploglik_known(spcov_initial_val, dispersion_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial_val)
    } else {
      cov_estimate_val <- use_laploglik(spcov_initial_val, dispersion_initial_val, data_object, estmethod,
                                        dist_matrix_list = dist_matrix_list, spcov_profiled = FALSE,
                                        randcov_initial = randcov_initial_val, randcov_profiled = FALSE,
                                        optim_dotlist = optim_dotlist
      )
    }
  }
}
