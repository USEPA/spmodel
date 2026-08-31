#' Get the (whitened) hat matrix (leverage)
#'
#' @param cov_betahat Covariance of fixed effects
#' @param cholprods A \code{cholprods} object
#'
#' @return Leverage (hat) values
#'
#' @noRd
get_hatvalues <- function(cov_betahat, SqrtSigInv_X) {
  # the hat matrix of the whitened residuals; only its diagonal (the
  # per-observation leverage) is needed, so get_diag_XVXt() forms it without
  # ever materializing the n x n product
  get_diag_XVXt(SqrtSigInv_X, cov_betahat)
}
