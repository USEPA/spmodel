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
  # Rather than refit the model n times with each observation dropped (which
  # would mean inverting a new (n-1)x(n-1) matrix each time), this uses a
  # partitioned-matrix (Sherman-Morrison-type) update: SigInv for the (n-1)
  # remaining observations can be recovered algebraically from the full-data
  # SigInv and the row/column for the dropped observation "obs".
  SigInv_mm <- SigInv[obs, obs] # a constant
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  newX <- Xmat[-obs, , drop = FALSE]
  newyX <- yX[-obs, , drop = FALSE]

  # SigInv %*% yX for the data with "obs" removed, obtained via the partitioned
  # inverse update rather than recomputing SigInv from scratch
  new_SigInv_oo_newyX <- SigInv_yX[-obs, , drop = FALSE] - SigInv_om %*% yX[obs, , drop = FALSE]
  newSigInv_newyX <- new_SigInv_oo_newyX - SigInv_om %*% (crossprod(SigInv_om, newyX) / SigInv_mm)

  newSigInv_newX <- newSigInv_newyX[, -1, drop = FALSE]
  newSigInv_newy <- newSigInv_newyX[, 1, drop = FALSE]
  # refit betahat using only the n-1 remaining observations
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, newSigInv_newX))))
  new_betahat <- new_covbetahat %*% crossprod(newX, newSigInv_newy)
  obs_c <- Sig[obs, -obs, drop = FALSE]
  # kriging predictor at "obs" using the refit betahat and the covariance
  # between "obs" and the remaining observations (universal kriging equation)
  new_pred <- Xmat[obs, , drop = FALSE] %*% new_betahat + obs_c %*% (newSigInv_newy - newSigInv_newX %*% new_betahat)

  # var
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c %*% newSigInv_newX
    new_SigInv_oo_obs_c <- tcrossprod(SigInv[-obs, -obs], obs_c) - SigInv_om %*% (crossprod(SigInv_om, t(obs_c)) / SigInv_mm)
    # kriging prediction variance: marginal variance minus variance explained
    # by the observed data, plus extra uncertainty from estimating betahat
    var_fit <- Sig[obs, obs] - obs_c %*% new_SigInv_oo_obs_c + Q %*% tcrossprod(new_covbetahat, Q)
    se_fit <- sqrt(var_fit)
  } else {
    se_fit <- NULL
  }

  # return
  list(pred = as.numeric(new_pred), se.fit = as.numeric(se_fit))
}

#' Get the exact (non-local) loocv standard error for iid errors
#'
#' @param obs An observation to leave out
#' @param cov_betahat The covariance matrix of the leave-one-out betahat
#' @param Xmat Model matrix
#' @param total_var The total (marginal) variance of a single observation
#'
#' @return A list with element \code{se.fit}, the loocv standard error
#'
#' @noRd
get_loocv_iid_se <- function(obs, cov_betahat, Xmat, total_var) {
  # for iid errors there's no spatial covariance to update via partitioned
  # inverse -- prediction variance is just the total variance plus the
  # variance of the fitted mean at this observation's covariates
  newX <- Xmat[obs, , drop = FALSE]
  var_fit <- newX %*% tcrossprod(cov_betahat, newX)
  se_fit <- sqrt(total_var + var_fit)
  list(se.fit = as.numeric(se_fit))
}
