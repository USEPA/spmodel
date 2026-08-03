#' Transform a dispersion parameter to the optimization (log) scale
#'
#' @param dispersion_initial A \code{dispersion_initial} object
#'
#' @return A list with \code{value}, the log-scale dispersion parameter
#'   (clamped to \eqn{[-50, 50]} for numerical stability), and \code{is_known},
#'   whether it is fixed rather than estimated
#'
#' @noRd
dispersion_orig2optim <- function(dispersion_initial) {
  # log transform moves the strictly-positive dispersion parameter onto an
  # unconstrained scale so optim() can search over all reals
  dispersion_orig2optim_val <- log(dispersion_initial$initial)
  names(dispersion_orig2optim_val) <- "dispersion_log"
  dispersion_orig2optim_is_known <- dispersion_initial$is_known
  names(dispersion_orig2optim_is_known) <- "dispersion_log"

  # return dispersion parameter vector
  # can consider lower limit like -10, 10 for numerical stability
  # clamp only values that are still being estimated -- known values are left
  # untouched since they are fixed by the user, not searched over
  dispersion_orig2optim_val <- ifelse(dispersion_orig2optim_val > 50 & !dispersion_orig2optim_is_known, 50, dispersion_orig2optim_val)
  dispersion_orig2optim_val <- ifelse(dispersion_orig2optim_val < -50 & !dispersion_orig2optim_is_known, -50, dispersion_orig2optim_val)

  dispersion_initial_list_val <- list(
    value = dispersion_orig2optim_val,
    is_known = dispersion_orig2optim_is_known
  )
  dispersion_initial_list_val
}
