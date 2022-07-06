#' Use semivariogram-weighted-least-squares for estimation
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param formula A formula
#' @param data Data
#' @param dist_matrix A distance matrix (Euclidean)
#' @param weights wls weights
#' @param optim_dotlist An optim dotlist
#' @param esv Empirical semivariogram
#'
#' @return The covariance parameter estimates
#'
#' @noRd
use_svloss <- function(spcov_initial, dist_matrix_list, esv, weights, optim_dotlist) {



  # transforming to optim paramters (log scale)
  spcov_orig2optim_val <- spcov_orig2optim(spcov_initial = spcov_initial, spcov_profiled = FALSE)

  # get optim par
  optim_par <- get_optim_par(spcov_orig2optim_val)

  # check optim dotlist
  optim_dotlist <- check_optim_method(optim_par, optim_dotlist)

  # performing optimization
  optim_output <- do.call("optim", c(
    list(
      par = optim_par,
      fn = svloss,
      spcov_orig2optim = spcov_orig2optim_val,
      esv = esv,
      weights = weights
    ),
    optim_dotlist
  ))

  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim_val, optim_output$par, spcov_profiled = FALSE)

  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim_val), spcov_orig_val = spcov_orig_val)

  # return parameter values and optim output
  optim_output <- list(
    method = optim_dotlist$method, control = optim_dotlist$control, value = optim_output$value,
    counts = optim_output$counts, convergence = optim_output$convergence,
    message = optim_output$message, hessian = if (optim_dotlist$hessian) optim_output$hessian else FALSE
  )
  # returning output
  list(
    spcov_params_val = spcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list, esv = esv,
    is_known = list(spcov = spcov_initial$is_known)
  )
}
