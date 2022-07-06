#' Get the (whitened) hat matrix (leverage)
#'
#' @param cov_betahat Covariance of fixed effects
#' @param cholprods A \code{cholprods} object
#'
#' @return Leverage (hat) values
#'
#' @noRd
get_hatvalues <- function(cov_betahat, SqrtSigInv_X) {
  # the hat matrix of the whitened residuals
  diag(SqrtSigInv_X %*% tcrossprod(cov_betahat, SqrtSigInv_X))
  # suggested cutoff is 2p / n or 3p / n (recall p / n is the average value) and sum
  # of hat matrix is p
}
