#' Draw exact conditional LM residuals with simulated coefficient uncertainty
#' @noRd
draw_conditional_exact_lm <- function(newdata, object, SqrtSigInv_y, SqrtSigInv_X, new_betahat, cov_lowchol_observed, samples, chunk_size = samples) {

  newdata_n <- NROW(newdata)

  cov_observed_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")

  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_observed, cov_observed_new)

  # ordinary kriging conditional covariance -- betahat uncertainty is already
  # supplied by the caller's simulated beta draws, so it is not added here
  cond_cov <- cov_new - crossprod(SqrtSigInv_c0, SqrtSigInv_c0)

  chol_cond_cov <- t(chol(cond_cov))

  # conditional mean: Sigma_{new,observed} Sigma_observed^-1 (y - X %*%
  # beta_b) == crossprod(SqrtSigInv_c0, SqrtSigInv_y) - crossprod(SqrtSigInv_c0,
  # SqrtSigInv_X) %*% beta_b
  cond_mu_y <- crossprod(SqrtSigInv_c0, SqrtSigInv_y)
  cond_mu_X <- crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  new_val <- matrix(NA_real_, newdata_n, samples)
  for (start in seq.int(1L, samples, by = chunk_size)) {
    cols <- seq.int(start, min(samples, start + chunk_size - 1L))
    cond_mu <- as.numeric(cond_mu_y) - cond_mu_X %*% new_betahat[, cols, drop = FALSE]
    new_val[, cols] <- chol_cond_cov %*% matrix(rnorm(newdata_n * length(cols)), newdata_n, length(cols)) + cond_mu
  }
  new_val
}

#' Draw exact conditional GLM residuals
#' @noRd
draw_conditional_exact_glm <- function(newdata, object, residual, cov_lowchol_observed, samples) {
  cov_observed_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_observed, cov_observed_new)
  wts_residual <- backsolve(t(cov_lowchol_observed), SqrtSigInv_c0)
  cond_cov <- as.matrix(cov_new - crossprod(SqrtSigInv_c0))
  chol_cond_cov <- t(chol(cond_cov))
  new_val <- chol_cond_cov %*% matrix(rnorm(NROW(newdata) * samples), NROW(newdata), samples)
  new_val + crossprod(wts_residual, residual)
}
