#' Get optim dotlist
#'
#' Get the dotlist for \code{optim()} within \code{splm()}
#'
#' @param ... Additional arguments to \code{optim()}
#'
#' @return An optim dotlist
#'
#' @noRd
get_optim_dotlist <- function(...) {
  # storing dotlist and setting defaults for optim
  # any optim() arguments the user passed via ... take precedence; only
  # fill in a default when the user did not already supply that argument
  dotlist <- list(...)

  # nelder-mead default with lower relative tolerance
  # derivative-free, so it works for arbitrary covariance functions without
  # requiring an analytic gradient of the (Laplace) log-likelihood
  if (!("method" %in% names(dotlist))) {
    dotlist$method <- "Nelder-Mead"
  }

  if (!("hessian" %in% names(dotlist))) {
    dotlist$hessian <- FALSE
  }

  if (!("control" %in% names(dotlist))) {
    dotlist$control <- list()
  }

  if (!("reltol" %in% names(dotlist$control))) {
    # prior to 0.13.0 reltol was 1e-4
    dotlist$control$reltol <- 1e-6
  }

  # lower/upper are unbounded here because covariance parameters are
  # optimized on a transformed ("optim") scale (see spcov_orig2optim /
  # dispersion_orig2optim) that maps constrained parameters (e.g. positive
  # variances) onto the whole real line, so optim() itself needs no bounds
  dotlist$lower <- -Inf
  dotlist$upper <- Inf

  # make optim dotlist
  # hardcode hessian false while developing satterthwaite
  dotlist$hessian <- FALSE
  optim_dotlist <- list(gr = NULL, method = dotlist$method, lower = dotlist$lower, upper = dotlist$upper, control = dotlist$control, hessian = dotlist$hessian)
}
