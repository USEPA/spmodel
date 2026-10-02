#' Prepare coefficient-first Gaussian simulation of latent and fixed effects
#'
#' @param object A fitted spglm or spgautor object.
#' @param cov_lowchol Optional observed covariance factor for an areal model.
#' @return Spatial and conditional-precision factors, coefficient
#'   covariance factor, coupling weights, and fitted centers in observation order.
#' @noRd
get_conditional_glm_joint <- function(object, cov_lowchol = NULL) {
  X <- model.matrix(object)
  eta <- fitted(object, type = "link")
  w <- w_offset_free(eta, model.offset(model.frame(object)))
  if (is.null(cov_lowchol)) {
    # Use the model covariance, including random effects and true partitions.
    # Local fitting groups supply estimates, not simulation covariance blocks.
    cov_lowchol <- t(chol(as.matrix(covmatrix(object))))
  }
  D <- get_D(object$family, eta, object$y, object$size,
    as.vector(coef(object, type = "dispersion")))
  SigInv <- chol2inv(t(cov_lowchol))
  SigInv_X <- SigInv %*% X
  diag(SigInv) <- diag(SigInv) - diag(D)
  cond_prec_upchol <- tryCatch(chol(SigInv), error = function(e) {
    stop("The conditional GLM latent precision is not positive definite; check the fitted model and its convergence.", call. = FALSE)
  })
  wts_latent <- backsolve(cond_prec_upchol,
    forwardsolve(t(cond_prec_upchol), SigInv_X))
  list(w = as.numeric(w), betahat = coef(object), wts_latent = wts_latent,
    residual_center = as.numeric(w - X %*% coef(object)), wts_residual = wts_latent - X,
    cov_lowchol = cov_lowchol, cond_prec_upchol = cond_prec_upchol,
    cov_betahat_lowchol = t(chol(vcov(object))))
}

#' Draw coefficients, then latent values conditional on those coefficients
#'
#' @param joint A context from get_conditional_glm_joint().
#' @param samples Number of draws.
#' @param residual Return latent residuals directly for spatial prediction.
#' @return Matrices w (or residual) and beta, with one draw per column.
#' @noRd
draw_conditional_glm_joint <- function(joint, samples, residual = FALSE) {
  delta <- joint$cov_betahat_lowchol %*%
    matrix(rnorm(length(joint$betahat) * samples), length(joint$betahat), samples)
  z <- backsolve(joint$cond_prec_upchol,
    matrix(rnorm(length(joint$w) * samples), length(joint$w), samples))
  beta <- sweep(delta, 1, joint$betahat, "+")
  rownames(beta) <- names(joint$betahat)
  if (residual) {
    return(list(residual = sweep(joint$wts_residual %*% delta + z,
      1, joint$residual_center, "+"), beta = beta))
  }
  list(w = sweep(joint$wts_latent %*% delta + z, 1, joint$w, "+"), beta = beta)
}
