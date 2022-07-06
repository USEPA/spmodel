#' Perform semivariogram weighted least squares estimation
#'
#' @param data_object The data object
#' @param formula A formula
#' @param spcov_initial The spatial initial object
#' @param estmethod Estimation method
#' @param weights sv-wls weights
#' @param optim_dotlist The optim dotlist
#' @param esv_dotlist The empirical semivariogram dotlist
#'
#' @return
#'
#' @noRd
cov_estimate_sv <- function(data_object, formula, spcov_initial, estmethod,
                            weights = weights, optim_dotlist, esv_dotlist) {

  # make NA spcov_initial
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy = data_object$anisotropy)

  # store distance matrix (if applicable)
  if (data_object$anisotropy) {
    dist_matrix_list <- NULL
  } else {
    dist_matrix_list <- lapply(data_object$obdata_list, function(x) spdist(x, data_object$xcoord, data_object$ycoord))
  }

  cov_initial_val <- cov_initial_search(spcov_initial_NA_val, estmethod, data_object,
    dist_matrix_list,
    weights = weights, esv_dotlist = esv_dotlist
  )

  spcov_initial_val <- cov_initial_val$spcov_initial_val

  if (data_object$anisotropy) {
    new_coords_list <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
      rotate = spcov_initial_val$initial[["rotate"]],
      scale = spcov_initial_val$initial[["scale"]]
    )
    dist_matrix_list <- lapply(new_coords_list, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))
  }


  if (all(spcov_initial_val$is_known)) {
    spcov_estimate_val <- use_svloss_known(spcov_initial_val, dist_matrix_list, cov_initial_val$esv, weights)
  } else {
    spcov_estimate_val <- use_svloss(spcov_initial_val, dist_matrix_list, cov_initial_val$esv, weights, optim_dotlist)
  }
  spcov_estimate_val
}
