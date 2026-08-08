#' Find relevant eigendecomposition-based quantities
#'
#' @param cov_matrix A covariance matrix
#' @param X A model matrix
#' @param y A response vector
#' @param ones A column vector of ones
#'
#' @return Relevant products of the covariance matrix's inverse and inverse
#'   square root (via its eigendecomposition) with \code{X}, \code{y}, and \code{ones};
#'   used instead of a Cholesky decomposition when the unique square root
#'   Sigma^{-1/2} must be computed (e.g., Pearson residuals in linear models).
#'
#' @noRd
get_eigenprods <- function(cov_matrix, X, y, ones) {
  # eigendecomposition Sigma = V diag(values) V' lets Sigma^{-1} and Sigma^{-1/2}
  # be formed by scaling the eigenvectors instead of factorizing Sigma directly
  eig <- eigen(Matrix::forceSymmetric(cov_matrix))
  # scaling each column of V by 1/eigenvalue (or 1/sqrt(eigenvalue)) is
  # equivalent to V %*% diag(1/values) but avoids building the diagonal matrix
  SigInv_p1 <- t(t(eig$vectors) * 1 / eig$values)
  SqrtSigInv_p1 <- t(t(eig$vectors) * 1 / sqrt(eig$values))
  # project X, y, ones into eigenvector space (V'X, V'y, V'ones); combined with
  # the scaled pieces above this reconstructs Sigma^{-1}X = V diag(1/values) V'X
  # without ever forming or inverting the full Sigma
  p2_X <- t(eig$vectors) %*% X
  p2_y <- t(eig$vectors) %*% y
  p2_ones <- t(eig$vectors) %*% ones
  SigInv_X <- SigInv_p1 %*% p2_X
  SigInv_y <- SigInv_p1 %*% p2_y
  SigInv_ones <- SigInv_p1 %*% p2_ones
  SqrtSigInv_X <- SqrtSigInv_p1 %*% p2_X
  SqrtSigInv_y <- SqrtSigInv_p1 %*% p2_y
  SqrtSigInv_ones <- SqrtSigInv_p1 %*% p2_ones
  list(
    SigInv_X = SigInv_X, SigInv_y = SigInv_y, SigInv_ones = SigInv_ones,
    SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_y = SqrtSigInv_y, SqrtSigInv_ones = SqrtSigInv_ones
  )
}

#' Parallel-friendly wrapper around \code{get_eigenprods()}
#'
#' @param cluster_list A list with elements \code{c} (covariance matrix),
#'   \code{x} (model matrix), \code{y} (response vector), and \code{o} (ones vector)
#'
#' @return The same value as \code{get_eigenprods()}, for use with \code{parallel::parLapply()}
#'
#' @noRd
get_eigenprods_parallel <- function(cluster_list) {
  cov_matrix <- cluster_list$c
  X <- cluster_list$x
  y <- cluster_list$y
  o <- cluster_list$o
  get_eigenprods(cov_matrix, X, y, o)
}
