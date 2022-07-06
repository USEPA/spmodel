#' Find relevant Cholesky quantities (upper, lower, and products)
#'
#' @param cov_matrix A covariance matrix
#' @param X A model matrix
#' @param y A response vector
#'
#' @return Relevant Cholesky quantities
#'
#' @noRd
get_cholprods <- function(cov_matrix, X, y) {
  Sig_upchol <- Matrix::chol(Matrix::forceSymmetric(cov_matrix))
  Sig_lowchol <- t(Sig_upchol)
  SqrtSigInv_X <- forwardsolve(Sig_lowchol, X)
  SqrtSigInv_y <- forwardsolve(Sig_lowchol, y)
  list(Sig_lowchol = Sig_lowchol, SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_y = SqrtSigInv_y)
}

get_cholprods_parallel <- function(cluster_list) {
  cov_matrix <- cluster_list$c
  X <- cluster_list$x
  y <- cluster_list$y
  get_cholprods(cov_matrix, X, y)
}
