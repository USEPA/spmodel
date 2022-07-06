#' Replace optim method with Brent if only one parameter requires optimization.
#'
#' @param optim_par A vector of parameters to optimize
#' @param optim_dotlist A dotlist with arguments for optim
#'
#' @return An edited method in \code{optim_dotlist} if \code{optim_par} has
#'   length 1. Brent method must have finite lower and upper bounds.
#'
#' @noRd
check_optim_method <- function(optim_par, optim_dotlist) {
  if (length(optim_par) == 1) {
    optim_dotlist$method <- "Brent"
    optim_dotlist$lower <- -50
    optim_dotlist$upper <- 50
  }
  optim_dotlist
}
