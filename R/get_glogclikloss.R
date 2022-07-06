#' Get composite (Gaussian log) likelhood loss
#'
#' @param spcov_params Spatial covariance parameters
#' @param residual_vector2 Squared residual vector
#' @param dist_vector Distance vector
#'
#' @return The composite (Gaussian log) likelhood loss
#'
#' @noRd
get_glogclikloss <- function(spcov_params, residual_vector2, dist_vector) {
  sigma2_val <- spcov_params[["de"]] + spcov_params[["ie"]]
  spcov_vec_val <- spcov_vector(spcov_params, dist_vector)
  sv_val <- sigma2_val - spcov_vec_val
  sum(residual_vector2 / (2 * sv_val) + log(sv_val))
}
