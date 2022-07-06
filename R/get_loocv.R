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
get_loocv <- function(obs, Sig, SigInv, Xmat, y, yX, SigInv_yX) {
  # SigInv_oo <- SigInv[-obs, -obs, drop = FALSE]
  SigInv_mm <- SigInv[obs, obs]
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  newX <- Xmat[-obs, , drop = FALSE]
  newyX <- yX[-obs, , drop = FALSE]

  # browser()
  new_SigInv_oo_newyX <- SigInv_yX[-obs, , drop = FALSE] - SigInv_om %*% yX[obs, , drop = FALSE]
  newSigInv_newyX <- new_SigInv_oo_newyX - SigInv_om %*% (crossprod(SigInv_om, newyX) / SigInv_mm)

  newSigInv_newX <- newSigInv_newyX[, -1, drop = FALSE]
  newSigInv_newy <- newSigInv_newyX[, 1, drop = FALSE]
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, newSigInv_newX))))
  new_betahat <- new_covbetahat %*% crossprod(newX, newSigInv_newy)
  obs_c <- Sig[obs, -obs]
  new_pred <- Xmat[obs, , drop = FALSE] %*% new_betahat + obs_c %*% (newSigInv_newy - newSigInv_newX %*% new_betahat)
  as.numeric(new_pred)

  # comp speed example
  # mx <- matrix(1:16, byrow = TRUE, nrow = 4)
  # mx2 <- matrix(8:1, byrow = TRUE, nrow = 4)
  ## remove 3
  # mx[-3, -3] %*% mx2[-3, ]
  # ((mx %*% mx2)[-3, ] - mx[-3, 3, drop = FALSE] %*% mx2[3, ])
  # SigInv_oo_newyX <- SigInv_oo_yX[-obs, ] - SigInv_oo[-obs, obs] %*% yX[obs, ]
  # this all avoids having to do the expensive, repeated SigInv[-obs, -obs, drop = FALSE] subset
  # and the expensive SigInv_oo %*% newyX computation
}



# 3-27-2022 update:
## this should work as above for local$method = "all" and work like prediction
## if local$method = "distance" or "covariance"
## take loocv out of model fitting process
## need inverse if local$method = "all"
## will likely need to return / use model.matrix.spmod and model.response.spmod
## and model.frame.spmod
