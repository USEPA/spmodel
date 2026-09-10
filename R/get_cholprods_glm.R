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
  # Cholesky factor Sigma = U'U (upper) where L = U'
  # (lower); using the triangular factor lets X and y be "whitened" via cheap
  # triangular solves instead of ever forming Sigma^{-1} explicitly
  Sig_upchol <- Matrix::chol(Matrix::forceSymmetric(cov_matrix))
  Sig_lowchol <- t(Sig_upchol)
  SqrtSigInv_X <- forwardsolve(Sig_lowchol, X)
  SqrtSigInv_y <- forwardsolve(Sig_lowchol, y)
  # unlike get_cholprods(), the GLM machinery (e.g. Laplace approximation)
  # needs the full inverse covariance and Sigma^{-1} %*% X explicitly, so
  # compute those here via the Cholesky factor rather than solve()
  SigInv <- chol2inv(Sig_upchol)
  SigInv_X <- backsolve(t(Sig_lowchol), SqrtSigInv_X)
  list(
    Sig_lowchol = Sig_lowchol, SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_y = SqrtSigInv_y,
    SigInv = SigInv, SigInv_X = SigInv_X
  )
}

#' Parallel-friendly wrapper around \code{get_cholprods_glm()}
#'
#' @param cluster_list A list with elements \code{c} (covariance matrix),
#'   \code{x} (model matrix), and \code{y} (response vector)
#'
#' @return The same value as \code{get_cholprods_glm()}, for use with \code{parallel::parLapply()}
#'
#' @noRd
get_cholprods_glm_parallel <- function(cluster_list) {
  cov_matrix <- cluster_list$c
  X <- cluster_list$x
  y <- cluster_list$y
  get_cholprods_glm(cov_matrix, X, y)
}
