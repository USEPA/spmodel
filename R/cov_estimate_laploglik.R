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
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy = data_object$anisotropy)

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

  } else {

  }

}
