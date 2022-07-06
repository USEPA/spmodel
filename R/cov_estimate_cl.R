#' Perform composite likelihood estimation for semivariogram marginal differences
#'   (Curriero and Lele, 1999)
#'
#'
#' @param data_object The data object
#' @param formula A formula
#' @param spcov_initial The spatial initial object
#' @param estmethod The estimation method
#' @param optim_dotlist The optim dotlist
#'
#' @return The composite log-likelihood estimates
#'
#' @noRd
#'
#' @references
#' Curriero, F. C., & Lele, S. (1999). A composite likelihood approach to
#'   semivariogram estimation. *Journal of Agricultural, biological, and
#'   Environmental statistics*, 9-28.
cov_estimate_cl <- function(data_object, formula, spcov_initial, estmethod, optim_dotlist) {


  # make NA spcov_initial
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy = data_object$anisotropy)

  # store distance matrix (if applicable)
  if (data_object$anisotropy) {
    dist_matrix_list <- NULL
  } else {
    dist_matrix_list <- lapply(data_object$obdata_list, function(x) spdist(x, data_object$xcoord, data_object$ycoord))
  }

  cov_initial_val <- cov_initial_search(
    spcov_initial_NA = spcov_initial_NA_val,
    estmethod = estmethod,
    data_object = data_object,
    dist_matrix_list = dist_matrix_list
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
    spcov_estimate_val <- use_glogclik_known(spcov_initial_val, data_object, dist_matrix_list, data_object$partition_list)
  } else {
    spcov_estimate_val <- use_glogclik(spcov_initial_val, data_object, dist_matrix_list, data_object$partition_list, optim_dotlist)
  }
  spcov_estimate_val
}
