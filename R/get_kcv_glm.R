#' Get the exact (non-local) kcv prediction and standard error for GLM-type models
#'
#' Block generalization of \code{\link{get_loocv_glm}()}: holds out a whole
#' fold (a vector of row indices) at once instead of a single observation, via
#' the same partitioned-matrix update with \code{solve()} against the
#' \code{m x m} blocks \code{SigInv[fold, fold]}/\code{mHinv[fold, fold]}
#' (\code{m} = fold size) in place of division by a scalar. At \code{m = 1}
#' this reduces to exactly \code{get_loocv_glm()}'s formula.
#'
#' @param fold A vector of row indices to leave out together
#' @param Sig The full covariance matrix
#' @param SigInv The full inverse covariance matrix
#' @param Xmat Model matrix
#' @param w The latent (link-scale) predictor vector
#' @param wX \code{cbind(w, Xmat)}
#' @param SigInv_wX \code{SigInv \%*\% wX}
#' @param mHinv The inverse of the negative Hessian of the Laplace log-likelihood
#' @param se.fit Whether to compute the standard error
#'
#' @return A list with elements \code{pred} (the link-scale kcv predictions
#'   for the fold) and \code{se.fit} (their standard errors, or \code{NULL} if
#'   \code{se.fit} is \code{FALSE}), computed via a partitioned-inverse
#'   (Sherman-Morrison-type) update rather than refitting the model with the
#'   fold removed
#'
#' @noRd
get_kcv_glm <- function(fold, Sig, SigInv, Xmat, w, wX, SigInv_wX, mHinv, se.fit) {
  SigInv_mm <- SigInv[fold, fold, drop = FALSE] # an m x m block (a scalar when m = 1)
  SigInv_om <- SigInv[-fold, fold, drop = FALSE]

  neww <- w[-fold, , drop = FALSE]
  newX <- Xmat[-fold, , drop = FALSE]

  # SigInv for the data with "fold" removed, via the partitioned inverse update
  new_SigInv <- SigInv[-fold, -fold] - SigInv_om %*% solve(SigInv_mm, t(SigInv_om))
  new_SigInv_newX <- new_SigInv %*% newX
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, new_SigInv_newX))))

  new_wts_beta <- tcrossprod(new_covbetahat, new_SigInv_newX)
  obs_c <- Sig[fold, -fold, drop = FALSE]
  obs_c_new_SigInv <- obs_c %*% new_SigInv
  obs_c_new_SigInv_newX <- obs_c %*% new_SigInv_newX
  # weights that map the remaining observations' latent w to the kriging
  # prediction at "fold" (universal kriging equation on the link scale)
  new_wts_pred <- Xmat[fold, , drop = FALSE] %*% new_wts_beta + obs_c_new_SigInv - obs_c_new_SigInv_newX %*% new_wts_beta

  new_pred <- new_wts_pred %*% neww

  # var
  if (se.fit) {
    Q <- Xmat[fold, , drop = FALSE] - obs_c_new_SigInv_newX
    var_fit <- Sig[fold, fold] - tcrossprod(obs_c_new_SigInv, obs_c) + Q %*% tcrossprod(new_covbetahat, Q)
    # mHinv (inverse negative Hessian of the Laplace loglik) captures the
    # extra uncertainty in w from the Laplace approximation itself; update it
    # the same way as SigInv above, then fold that uncertainty into var_fit
    mHinv_mm <- mHinv[fold, fold, drop = FALSE]
    mHinv_om <- mHinv[-fold, fold, drop = FALSE]
    newmHinv <- mHinv[-fold, -fold] - mHinv_om %*% solve(mHinv_mm, t(mHinv_om))
    var_adj <- var_fit + new_wts_pred %*% tcrossprod(newmHinv, new_wts_pred)
    # only the diagonal (marginal SE per point) is returned, matching
    # get_loocv_glm()'s per-observation output shape
    se_fit <- sqrt(diag(var_adj))
  } else {
    se_fit <- NULL
  }

  # return
  list(pred = as.numeric(new_pred), se.fit = se_fit)
}
