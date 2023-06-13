get_loocv_glm <- function(obs, Sig, SigInv, Xmat, w, wX, SigInv_wX, mHinv, se.fit) {
  SigInv_mm <- SigInv[obs, obs] # a constant
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  neww <- w[-obs, , drop = FALSE]
  newX <- Xmat[-obs, , drop = FALSE]
  newwX <- wX[-obs, , drop = FALSE]

  new_SigInv <- SigInv[-obs, -obs] - tcrossprod(SigInv_om, SigInv_om) / SigInv_mm
  new_SigInv_newX <- new_SigInv %*% newX
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, new_SigInv_newX))))

  new_wts_beta <- tcrossprod(new_covbetahat, new_SigInv_newX)
  obs_c <- Sig[obs, -obs, drop = FALSE]
  obs_c_new_SigInv <- obs_c %*% new_SigInv
  obs_c_new_SigInv_newX <- obs_c %*% new_SigInv_newX
  new_wts_pred <- Xmat[obs, , drop = FALSE] %*% new_wts_beta + obs_c %*% new_SigInv - obs_c_new_SigInv_newX %*% new_wts_beta


  new_pred <- new_wts_pred %*% neww

  # var
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c_new_SigInv_newX
    var_fit <- Sig[obs, obs] - tcrossprod(obs_c_new_SigInv, obs_c) + Q %*% tcrossprod(new_covbetahat, Q)
    mHinv_mm <- mHinv[obs, obs]
    mHinv_om <- mHinv[-obs, obs, drop = FALSE]
    newmHinv <- mHinv[-obs, -obs] - tcrossprod(mHinv_om, mHinv_om) / mHinv_mm
    var_adj <- as.numeric(var_fit + new_wts_pred %*% tcrossprod(newmHinv, new_wts_pred))
    se_fit <- sqrt(var_adj)
  } else {
    se_fit <- NULL
  }

  # return
  list(pred = as.numeric(new_pred), se.fit = as.numeric(se_fit))
}
