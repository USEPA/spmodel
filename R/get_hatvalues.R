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
  # only the diagonal (per-observation leverage) is needed, so the full
  # n x n hat matrix is never formed -- diag(A %*% B %*% t(A)) is computed via
  # tcrossprod/diag rather than materializing the dense product
  diag(SqrtSigInv_X %*% tcrossprod(cov_betahat, SqrtSigInv_X))
}
