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
  dotlist <- list(...)

  ## l-bfgs-b deafult
  # if (!("method" %in% names(dotlist))) {
  #   dotlist$method <- "L-BFGS-B"
  # }

  # nelder-mead default with lower relative tolerance
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
    dotlist$control$reltol <- 1e-4
  }

  dotlist$lower <- -Inf
  dotlist$upper <- Inf

  # make optim dotlist
  optim_dotlist <- list(gr = NULL, method = dotlist$method, lower = dotlist$lower, upper = dotlist$upper, control = dotlist$control, hessian = dotlist$hessian)
}
