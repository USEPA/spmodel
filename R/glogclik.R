#' Find (Gaussian log) composite likelihood while optimzing
#'
#' @param par parameters to optimize over
#' @param spcov_orig2optim a \code{spcov_orig2optim} object
#' @param residual_vector2 Squared residuals
#' @param dist_vector Distance vector
#'
#' @return (Gaussian log) composite likelihood while optimzing
#'
#' @noRd
glogclik <- function(par, spcov_orig2optim, residual_vector2, dist_vector, data_object) {
  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim, par, spcov_profiled = FALSE,
                                     data_object = data_object)
  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim), spcov_orig_val = spcov_orig_val)
  glogclikloss_val <- get_glogclikloss(spcov_params_val, residual_vector2, dist_vector)
}
