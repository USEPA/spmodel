#' Use semivariogram-weighted-least-squares for estimation with known covariance parameters
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param formula A formula
#' @param data Data
#' @param dist_matrix A distance matrix (Euclidean)
#' @param weights wls weights
#' @param esv Empirical semivariogram
#'
#' @return The known covariance parameters
#'
#' @noRd
use_svloss_known <- function(spcov_initial, dist_matrix_list, esv, weights) {


  # get covariance parameters
  spcov_params_val <- get_spcov_params(class(spcov_initial), spcov_initial$initial)
  # get sv loss
  svloss_val <- get_svloss(
    spcov_params_val,
    esv = esv,
    weights = weights
  )

  # return parameter values and optim output
  optim_output <- list(
    method = NA, control = NA, value = svloss_val,
    counts = NA, convergence = NA,
    message = NA, hessian = NA
  )
  # returning output
  list(
    spcov_params_val = spcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list, esv = esv,
    is_known = list(spcov = spcov_initial$is_known)
  )
}
