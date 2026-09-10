#' Find relevant eigendecomposition-based quantities (with the full inverse)
#'
#' @param cov_matrix A covariance matrix
#' @param X A model matrix
#' @param y A response vector
#' @param ones A column vector of ones
#'
#' @return The same quantities as \code{get_eigenprods()}, plus \code{SigInv}, the
#'   full covariance matrix inverse
#'
#' @noRd
get_eigenprods_glm <- function(cov_matrix, X, y, ones) {
  # same eigendecomposition derivation as get_eigenprods()
  eig <- eigen(Matrix::forceSymmetric(cov_matrix))
  SigInv_p1 <- t(t(eig$vectors) * 1 / eig$values)
  SqrtSigInv_p1 <- t(t(eig$vectors) * 1 / sqrt(eig$values))
  p2_X <- t(eig$vectors) %*% X
  p2_y <- t(eig$vectors) %*% y
  p2_ones <- t(eig$vectors) %*% ones
  SigInv_X <- SigInv_p1 %*% p2_X
  SigInv_y <- SigInv_p1 %*% p2_y
  SigInv_ones <- SigInv_p1 %*% p2_ones
  SqrtSigInv_X <- SqrtSigInv_p1 %*% p2_X
  SqrtSigInv_y <- SqrtSigInv_p1 %*% p2_y
  SqrtSigInv_ones <- SqrtSigInv_p1 %*% p2_ones
  # unlike get_eigenprods(), the full inverse Sigma^{-1} = V diag(1/values) V' is
  # also formed here (needed downstream for GLM working-response BLUP updates)
  SigInv <- SigInv_p1 %*% t(eig$vectors)
  list(
    SigInv_X = SigInv_X, SigInv_y = SigInv_y, SigInv_ones = SigInv_ones,
    SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_y = SqrtSigInv_y, SqrtSigInv_ones = SqrtSigInv_ones,
    SigInv = SigInv
  )
}

#' Parallel-friendly wrapper around \code{get_eigenprods_glm()}
#'
#' @param cluster_list A list with elements \code{c} (covariance matrix),
#'   \code{x} (model matrix), \code{y} (response vector), and \code{o} (ones vector)
#'
#' @return The same value as \code{get_eigenprods_glm()}, for use with \code{parallel::parLapply()}
#'
#' @noRd
get_eigenprods_glm_parallel <- function(cluster_list) {
  cov_matrix <- cluster_list$c
  X <- cluster_list$x
  y <- cluster_list$y
  o <- cluster_list$o
  get_eigenprods_glm(cov_matrix, X, y, o)
}
