#' Warn when the covariance-parameter optimizer did not converge
#'
#' \code{optim()}'s reported convergence code is computed for every
#' optimization-based fit but was previously never inspected. \code{NA} means
#' \code{optim()} was never called (all covariance parameters were fixed/known,
#' or a closed-form shortcut applies), which is not a failure and must not warn.
#'
#' @noRd
warn_optim_convergence <- function(convergence) {
  if (is.na(convergence) || convergence == 0) {
    return(invisible())
  }
  warning(
    paste0(
      "optim() did not converge while fitting the model ",
      "(convergence code ", convergence, "). Fitted model estimates may be unreliable."
    ),
    call. = FALSE
  )
  invisible()
}