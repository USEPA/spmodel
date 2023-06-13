#' Find relevant Cholesky quantities (upper, lower, products, and inverse)
#'
#' @param cov_matrix A covariance matrix
#' @param X A model matrix
#' @param y A response vector
#'
#' @return Relevant Cholesky quantities
#'
#' @noRd
get_cholprods_glm <- function(cov_matrix, X, y) {
  Sig_upchol <- Matrix::chol(Matrix::forceSymmetric(cov_matrix))
  Sig_lowchol <- t(Sig_upchol)
  SqrtSigInv_X <- forwardsolve(Sig_lowchol, X)
  SqrtSigInv_y <- forwardsolve(Sig_lowchol, y)
  SigInv <- chol2inv(Sig_upchol)
  SigInv_X <- backsolve(t(Sig_lowchol), SqrtSigInv_X)
  # list(Sig_lowchol = Sig_lowchol, SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_y = SqrtSigInv_y)
  list(
    Sig_lowchol = Sig_lowchol, SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_y = SqrtSigInv_y,
    SigInv = SigInv, SigInv_X = SigInv_X
  )
}

get_cholprods_glm_parallel <- function(cluster_list) {
  cov_matrix <- cluster_list$c
  X <- cluster_list$x
  y <- cluster_list$y
  get_cholprods_glm(cov_matrix, X, y)
}
