#' Compute sv-wls loss
#'
#' @param par Parameter vector to optimize over
#' @param spcov_orig2optim A \code{spcov_opig2optim} object
#' @param gamma Empirical semivariogram value
#' @param weights Empirical semivariogram weights
#' @param dist_vector Distance vector
#' @param np Empirical semivariogram distance bins
#'
#' @return svlw-s loss
#'
#' @noRd
svloss <- function(par, spcov_orig2optim, esv, weights, data_object) {
  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim, par, spcov_profiled = FALSE,
                                     data_object = data_object)
  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim), spcov_orig_val = spcov_orig_val)
  svloss_val <- get_svloss(spcov_params_val, esv, weights)
}
