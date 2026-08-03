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
  # sigma2_val is the total variance (sill): dependent + independent error
  sigma2_val <- spcov_params[["de"]] + spcov_params[["ie"]]
  # covariance at each pairwise distance, under the candidate spcov_params
  spcov_vec_val <- spcov_vector(spcov_params, dist_vector)
  # sv_val is the (semi)variance implied by the candidate parameters at each
  # distance -- variance of a pairwise difference is sill minus covariance
  sv_val <- sigma2_val - spcov_vec_val
  # weighted least-squares-style loss: penalizes both the squared mismatch
  # between observed and model-implied variance (scaled by model variance) and
  # the model variance itself, so candidate parameters are scored like a
  # composite (pairwise) Gaussian log-likelihood without needing a full n x n
  # covariance matrix
  sum(residual_vector2 / (2 * sv_val) + log(sv_val))
}
