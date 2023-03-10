get_eigenprods_glm <- function(cov_matrix, X, y, ones) {

  eig <- eigen(Matrix::forceSymmetric(cov_matrix))
  SigInv_p1 <- t(t(eig$vectors) * 1/eig$values)
  SqrtSigInv_p1 <- t(t(eig$vectors) * 1/sqrt(eig$values))
  p2_X <- t(eig$vectors) %*% X
  p2_y <- t(eig$vectors) %*% y
  p2_ones <- t(eig$vectors) %*% ones
  SigInv_X <- SigInv_p1 %*% p2_X
  SigInv_y <- SigInv_p1 %*% p2_y
  SigInv_ones <- SigInv_p1 %*% p2_ones
  SqrtSigInv_X <- SqrtSigInv_p1 %*% p2_X
  SqrtSigInv_y <- SqrtSigInv_p1 %*% p2_y
  SqrtSigInv_ones <- SqrtSigInv_p1 %*% p2_ones
  SigInv <- SigInv_p1 %*% t(eig$vectors)
  list(SigInv_X = SigInv_X, SigInv_y = SigInv_y, SigInv_ones = SigInv_ones,
       SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_y = SqrtSigInv_y, SqrtSigInv_ones = SqrtSigInv_ones,
       SigInv = SigInv)
}

get_eigenprods_glm_parallel <- function(cluster_list) {
  cov_matrix <- cluster_list$c
  X <- cluster_list$x
  y <- cluster_list$y
  o <- cluster_list$o
  get_eigenprods(cov_matrix, X, y, o)
}
