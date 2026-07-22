get_conditional_new_from_base <- function(newdata, object, base_val, cov_lowchol_base, samples) {

  newdata_n <- NROW(newdata)
  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")


  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, cov_base_new)
  SqrtSigInv_base_val <- forwardsolve(cov_lowchol_base, base_val)

  cond_cov <- cov_new - crossprod(SqrtSigInv_c0, SqrtSigInv_c0)
  chol_cond_cov <- t(chol(cond_cov))
  new_val <- vapply(seq_len(samples), function(x) as.numeric(chol_cond_cov %*% rnorm(newdata_n)), numeric(newdata_n))

  cond_mu <- crossprod(SqrtSigInv_c0, SqrtSigInv_base_val)
  new_val <- new_val + cond_mu
}

get_conditional_new_from_base_adjust <- function(newdata_list, object, base_val, cov_lowchol_base, samples, SqrtSigInv_X, cov_betahat) {

  x0 <- newdata_list$x0
  newdata <- newdata_list$newdata

  newdata_n <- NROW(newdata)

  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")


  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, cov_base_new)
  SqrtSigInv_base_val <- forwardsolve(cov_lowchol_base, base_val)

  cond_cov <- cov_new - crossprod(SqrtSigInv_c0, SqrtSigInv_c0)
  H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  cond_cov <- cond_cov + H %*% Matrix::tcrossprod(cov_betahat, H)

  spcov_val <- coef(object, type = "spcov")
  if (spcov_val[["de"]] == 0 && is.null(coef(object, type = "randcov"))) {
    chol_cond_cov <- Matrix::Diagonal(NROW(cond_cov))
    diag(chol_cond_cov) <- sqrt(diag(chol_cond_cov))
  } else {
    chol_cond_cov <- t(chol(cond_cov))
  }

  new_val <- vapply(seq_len(samples), function(x) as.numeric(chol_cond_cov %*% rnorm(newdata_n)), numeric(newdata_n))

  cond_mu <- crossprod(SqrtSigInv_c0, SqrtSigInv_base_val)
  new_val <- new_val + cond_mu
}

get_conditional_new_from_base_adjust_glm <- function(newdata_list, object, base_val, cov_lowchol_base, samples, SqrtSigInv_X, cov_betahat, SigInv, SigInv_X, wts_beta, cov_lowchol_mH) {

  x0 <- newdata_list$x0
  newdata <- newdata_list$newdata

  newdata_n <- NROW(newdata)
  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")


  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, cov_base_new)
  SqrtSigInv_base_val <- forwardsolve(cov_lowchol_base, base_val)

  cond_cov <- cov_new - crossprod(SqrtSigInv_c0, SqrtSigInv_c0)
  H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  cond_cov <- cond_cov + H %*% Matrix::tcrossprod(cov_betahat, H)

  c0 <- t(cov_base_new)
  wts_pred <- x0 %*% wts_beta + c0 %*% SigInv - (c0 %*% SigInv_X) %*% wts_beta
  wts_pred <- t(wts_pred)
  SqrtmHInv_wts_pred <- forwardsolve(cov_lowchol_mH, wts_pred)
  var_adj <- crossprod(SqrtmHInv_wts_pred, SqrtmHInv_wts_pred)


  cov_cov <- var_adj + cond_cov

  chol_cond_cov <- t(chol(cond_cov))
  new_val <- vapply(seq_len(samples), function(x) as.numeric(chol_cond_cov %*% rnorm(newdata_n)), numeric(newdata_n))

  cond_mu <- crossprod(SqrtSigInv_c0, SqrtSigInv_base_val)
  new_val <- new_val + cond_mu
}
