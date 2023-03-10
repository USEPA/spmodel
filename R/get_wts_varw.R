get_wts_varw <- function(family, Xmat, y, w, size, dispersion, cov_lowchol, x0, c0) {

  SigInv <- chol2inv(t(cov_lowchol)) # works on upchol
  SigInv_X <- SigInv %*% Xmat
  cov_betahat <- chol2inv(chol(crossprod(Xmat, SigInv_X)))
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  Ptheta <- SigInv - SigInv_X %*% wts_beta

  d <-  get_d(family, w, y, size, dispersion)
  # and then the gradient vector
  # g <-  d - Ptheta %*% w
  # Next, compute H
  D <- get_D(family, w, y, size, dispersion)
  H <- D - Ptheta
  mHInv <- solve(-H) #chol2inv(chol(Matrix::forceSymmetric(-H))) # solve(-H)

  if (NROW(x0) == 1) {
    wts_pred <- x0 %*% wts_beta + c0 %*% SigInv - (c0 %*% SigInv_X) %*% wts_beta
    var_adj <- as.numeric(wts_pred %*% tcrossprod(mHInv, wts_pred))
  } else {
    var_adj <- vapply(seq_len(NROW(x0)), function(x) {
      x0_new <- x0[x, , drop = FALSE]
      c0_new <- c0[x, , drop = FALSE]
      wts_pred <- x0_new %*% wts_beta + c0_new %*% SigInv - (c0_new %*% SigInv_X) %*% wts_beta
      as.numeric(wts_pred %*% tcrossprod(mHInv, wts_pred))
    }, numeric(1))
  }
  var_adj
}
