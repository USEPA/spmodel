#' Get the exact (non-local) loocv prediction and standard error for GLM-type models
#'
#' @param obs An observation to leave out
#' @param Sig The full covariance matrix
#' @param SigInv The full inverse covariance matrix
#' @param Xmat Model matrix
#' @param w The latent (link-scale) predictor vector
#' @param wX \code{cbind(w, Xmat)}
#' @param SigInv_wX \code{SigInv \%*\% wX}
#' @param mHinv The inverse of the negative Hessian of the Laplace log-likelihood
#' @param se.fit Whether to compute the standard error
#'
#' @return A list with elements \code{pred} (the link-scale loocv prediction)
#'   and \code{se.fit} (its standard error, or \code{NULL} if \code{se.fit} is \code{FALSE}),
#'   computed via a partitioned-inverse (Sherman-Morrison-type) update rather
#'   than refitting the model with the observation removed
#'
#' @noRd
get_loocv_glm <- function(obs, Sig, SigInv, Xmat, w, wX, SigInv_wX, mHinv, se.fit) {
  # like get_loocv(), this avoids literally refitting with "obs" dropped by
  # updating SigInv via a partitioned-matrix (Sherman-Morrison-type) formula;
  # the GLM case works on the link-scale latent predictor w (from the Laplace
  # approximation) rather than the response y directly
  SigInv_mm <- SigInv[obs, obs] # a constant
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  neww <- w[-obs, , drop = FALSE]
  newX <- Xmat[-obs, , drop = FALSE]
  newwX <- wX[-obs, , drop = FALSE]

  # SigInv for the data with "obs" removed, via the partitioned inverse update
  new_SigInv <- SigInv[-obs, -obs] - tcrossprod(SigInv_om, SigInv_om) / SigInv_mm
  new_SigInv_newX <- new_SigInv %*% newX
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, new_SigInv_newX))))

  new_wts_beta <- tcrossprod(new_covbetahat, new_SigInv_newX)
  obs_c <- Sig[obs, -obs, drop = FALSE]
  obs_c_new_SigInv <- obs_c %*% new_SigInv
  obs_c_new_SigInv_newX <- obs_c %*% new_SigInv_newX
  # weights that map the remaining observations' latent w to the kriging
  # prediction at "obs" (universal kriging equation on the link scale)
  new_wts_pred <- Xmat[obs, , drop = FALSE] %*% new_wts_beta + obs_c %*% new_SigInv - obs_c_new_SigInv_newX %*% new_wts_beta


  new_pred <- new_wts_pred %*% neww

  # var
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c_new_SigInv_newX
    var_fit <- Sig[obs, obs] - tcrossprod(obs_c_new_SigInv, obs_c) + Q %*% tcrossprod(new_covbetahat, Q)
    # mHinv (inverse negative Hessian of the Laplace loglik) captures the
    # extra uncertainty in w from the Laplace approximation itself; update it
    # the same way as SigInv above, then fold that uncertainty into var_fit
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
