#' Prepare the fitted joint Gaussian approximation for latent and fixed effects
#'
#' @param object A fitted spglm object.
#' @return Lower Cholesky factors of the negative latent Hessian and coefficient
#'   covariance, coupling weights, and fitted centers in original observation order.
#' @noRd
get_conditional_glm_joint <- function(object) {
  X <- model.matrix(object)
  eta <- fitted(object, type = "link")
  w <- w_offset_free(eta, model.offset(model.frame(object)))
  index <- object$local_index
  if (is.null(index) && !is.null(object$partition_factor)) {
    index <- model.frame(object$partition_factor, object$obdata)[[1L]]
  }
  if (is.null(index)) index <- rep(1L, NROW(X))
  groups <- split(seq_len(NROW(X)), index, drop = TRUE)
  SigInv <- matrix(0, NROW(X), NROW(X))
  for (rows in groups) {
    block <- object
    block$obdata <- object$obdata[rows, , drop = FALSE]
    SigInv[rows, rows] <- chol2inv(chol(as.matrix(covmatrix(block))))
  }
  SigInv_X <- SigInv %*% X
  invcov_betahat_sum <- crossprod(X, SigInv_X)
  # Match the regularized inverse used by fitting's multi-block SMW Hessian.
  if (length(groups) > 1L) diag(invcov_betahat_sum) <- diag(invcov_betahat_sum) + object$diagtol
  cov_betahat_noadjust <- chol2inv(chol(invcov_betahat_sum))
  cov_betahat <- vcov(object, var_correct = FALSE)
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  Ptheta <- as.matrix(SigInv - SigInv_X %*% tcrossprod(cov_betahat_noadjust, SigInv_X))
  D <- get_D(object$family, eta, object$y, object$size,
    as.vector(coef(object, type = "dispersion")))
  H <- D - Ptheta
  cov_lowchol_mH <- tryCatch(t(chol(as.matrix(-H))), error = function(e) {
    stop("The conditional GLM latent precision is not positive definite; check the fitted model and its convergence.", call. = FALSE)
  })
  list(w = as.numeric(w), betahat = coef(object), wts_beta = wts_beta,
    cov_lowchol_mH = cov_lowchol_mH, cov_betahat_lowchol = t(chol(cov_betahat)))
}

#' Draw coupled latent and fixed effects
#'
#' @param joint A context from get_conditional_glm_joint().
#' @param samples Number of draws.
#' @return Matrices w (offset-free) and beta, with one draw per column.
#' @noRd
draw_conditional_glm_joint <- function(joint, samples) {
  u <- backsolve(t(joint$cov_lowchol_mH),
    matrix(rnorm(length(joint$w) * samples), length(joint$w), samples))
  new_betahat <- joint$wts_beta %*% u + joint$cov_betahat_lowchol %*%
    matrix(rnorm(length(joint$betahat) * samples), length(joint$betahat), samples)
  new_betahat <- sweep(new_betahat, 1, joint$betahat, "+")
  rownames(new_betahat) <- names(joint$betahat)
  list(w = sweep(u, 1, joint$w, "+"), beta = new_betahat)
}
