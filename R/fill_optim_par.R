#' Fill parameters being optimized with known parameters
#'
#' @param cov_orig2optim A \code{cov_orig2optim} object (parameters on optim scale)
#' @param par The parameters to optimize over in optim
#'
#' @return A covariance parameter vector (on the optim scale)
#'
#' @noRd
fill_optim_par <- function(cov_orig2optim, par) {
  cov_orig2optim$value[!cov_orig2optim$is_known] <- par
  cov_orig2optim$value
}
