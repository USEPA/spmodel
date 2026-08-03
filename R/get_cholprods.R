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
  # Cholesky factor Sigma = U'U (upper) where L = U'
  # (lower); using the triangular factor lets X and y be "whitened" via cheap
  # triangular solves instead of ever forming Sigma^{-1} explicitly
  Sig_upchol <- Matrix::chol(Matrix::forceSymmetric(cov_matrix))
  Sig_lowchol <- t(Sig_upchol)
  # forwardsolve(L, b) computes L^{-1} b == a by La = b and
  # solving for a (generalized least
  # squares whitening), reused throughout likelihood and prediction code
  # then (L^{-1} b)^t (L^{-1} b) = b^t Sigma^{-1} b = a^t a
  SqrtSigInv_X <- forwardsolve(Sig_lowchol, X)
  SqrtSigInv_y <- forwardsolve(Sig_lowchol, y)
  list(Sig_lowchol = Sig_lowchol, SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_y = SqrtSigInv_y)
}

#' Parallel-friendly wrapper around \code{get_cholprods()}
#'
#' @param cluster_list A list with elements \code{c} (covariance matrix),
#'   \code{x} (model matrix), and \code{y} (response vector)
#'
#' @return The same value as \code{get_cholprods()}, for use with \code{parallel::parLapply()}
#'
#' @noRd
get_cholprods_parallel <- function(cluster_list) {
  cov_matrix <- cluster_list$c
  X <- cluster_list$x
  y <- cluster_list$y
  get_cholprods(cov_matrix, X, y)
}
