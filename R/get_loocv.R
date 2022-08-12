#' Get loocv residual
#'
#' @param obs An observation to leave out
#' @param Sig The full covariance matrix
#' @param SigInv The full inverse covariance matrix
#' @param Xmat Model matrix
#' @param y response vector
#'
#' @return A loocv residual
#'
#' @noRd
get_loocv <- function(obs, Sig, SigInv, Xmat, y, yX, SigInv_yX, se.fit) {
  SigInv_mm <- SigInv[obs, obs] # a constant
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  newX <- Xmat[-obs, , drop = FALSE]
  newyX <- yX[-obs, , drop = FALSE]

  new_SigInv_oo_newyX <- SigInv_yX[-obs, , drop = FALSE] - SigInv_om %*% yX[obs, , drop = FALSE]
  newSigInv_newyX <- new_SigInv_oo_newyX - SigInv_om %*% (crossprod(SigInv_om, newyX) / SigInv_mm)

  newSigInv_newX <- newSigInv_newyX[, -1, drop = FALSE]
  newSigInv_newy <- newSigInv_newyX[, 1, drop = FALSE]
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, newSigInv_newX))))
  new_betahat <- new_covbetahat %*% crossprod(newX, newSigInv_newy)
  obs_c <- Sig[obs, -obs, drop = FALSE]
  new_pred <- Xmat[obs, , drop = FALSE] %*% new_betahat + obs_c %*% (newSigInv_newy - newSigInv_newX %*% new_betahat)

  # var
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c %*% newSigInv_newX
    new_SigInv_oo_obs_c <- tcrossprod(SigInv[-obs, -obs], obs_c) - SigInv_om %*% (crossprod(SigInv_om, t(obs_c)) / SigInv_mm)
    var_fit <- Sig[obs, obs] - obs_c %*% new_SigInv_oo_obs_c + Q %*% tcrossprod(new_covbetahat, Q)
    se_fit <- sqrt(var_fit)
  } else {
    se_fit <- NULL
  }

  # return
  list(pred = as.numeric(new_pred), se.fit = as.numeric(se_fit))
}
