#' Get kcv fold residual
#'
#' Block generalization of \code{\link{get_loocv}()}: rather than holding out a
#' single observation \code{obs}, holds out a whole fold (a vector of row
#' indices) at once via the same partitioned-matrix (Sherman-Morrison-type)
#' update, replacing division by the scalar \code{SigInv[obs, obs]} with
#' \code{solve()} against the \code{m x m} block \code{SigInv[fold, fold]}
#' (\code{m} = fold size). At \code{m = 1} this reduces to exactly
#' \code{get_loocv()}'s formula.
#'
#' @param fold A vector of row indices to leave out together
#' @param Sig The full covariance matrix
#' @param SigInv The full inverse covariance matrix
#' @param Xmat Model matrix
#' @param y response vector
#'
#' @return A kcv residual
#'
#' @noRd
get_kcv <- function(fold, Sig, SigInv, Xmat, y, yX, SigInv_yX, se.fit) {
  SigInv_mm <- SigInv[fold, fold, drop = FALSE] # an m x m block (a scalar when m = 1)
  SigInv_om <- SigInv[-fold, fold, drop = FALSE]

  newX <- Xmat[-fold, , drop = FALSE]
  newyX <- yX[-fold, , drop = FALSE]

  # SigInv %*% yX for the data with "fold" removed, obtained via the partitioned
  # inverse update rather than recomputing SigInv from scratch
  new_SigInv_oo_newyX <- SigInv_yX[-fold, , drop = FALSE] - SigInv_om %*% yX[fold, , drop = FALSE]
  newSigInv_newyX <- new_SigInv_oo_newyX - SigInv_om %*% solve(SigInv_mm, crossprod(SigInv_om, newyX))

  newSigInv_newX <- newSigInv_newyX[, -1, drop = FALSE]
  newSigInv_newy <- newSigInv_newyX[, 1, drop = FALSE]
  # refit betahat using only the remaining observations
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, newSigInv_newX))))
  new_betahat <- new_covbetahat %*% crossprod(newX, newSigInv_newy)
  obs_c <- Sig[fold, -fold, drop = FALSE]
  # kriging predictor at "fold" using the refit betahat and the covariance
  # between the fold and the remaining observations (universal kriging equation)
  new_pred <- Xmat[fold, , drop = FALSE] %*% new_betahat + obs_c %*% (newSigInv_newy - newSigInv_newX %*% new_betahat)

  # var
  if (se.fit) {
    Q <- Xmat[fold, , drop = FALSE] - obs_c %*% newSigInv_newX
    new_SigInv_oo_obs_c <- tcrossprod(SigInv[-fold, -fold], obs_c) - SigInv_om %*% solve(SigInv_mm, crossprod(SigInv_om, t(obs_c)))
    # kriging prediction covariance among the fold's own points: marginal
    # covariance minus variance explained by the observed data, plus extra
    # uncertainty from estimating betahat -- only the diagonal (marginal SE per
    # point) is returned, matching get_loocv()'s per-observation output shape
    var_fit <- Sig[fold, fold] - obs_c %*% new_SigInv_oo_obs_c + Q %*% tcrossprod(new_covbetahat, Q)
    se_fit <- sqrt(diag(var_fit))
  } else {
    se_fit <- NULL
  }

  # return
  list(pred = as.numeric(new_pred), se.fit = se_fit)
}
