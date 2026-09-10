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
  # mark which spatial covariance parameters the user fixed (vs. left to be
  # estimated) so the grid search below only searches over the free ones
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy = data_object$anisotropy)

  # store distance matrix (if applicable); skipped under anisotropy since distances
  # depend on the rotate/scale parameters, which aren't known until estimation runs
  if (data_object$anisotropy) {
    dist_matrix_list <- NULL
  } else {
    dist_matrix_list <- lapply(data_object$obdata_list, function(x) spdist(x, data_object$xcoord, data_object$ycoord))
  }

  # grid search for good optimizer starting values (avoids composite-likelihood
  # optimization landing in a poor local optimum)
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

  # choose known-evaluation vs. optimization -- see run_cl_dispatch() in
  # cov_estimate_dispatch_helpers.R
  spcov_estimate_val <- run_cl_dispatch(spcov_initial_val, data_object, dist_matrix_list, optim_dotlist)
  spcov_estimate_val
}
