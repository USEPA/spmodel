#' Compute the prediction variance adjustment for GLM-type models
#'
#' @param family The response family
#' @param Xmat Model matrix
#' @param y Response vector (unused directly; \code{w} carries the fitted latent values)
#' @param w The latent (link-scale) predictor vector
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter
#' @param cov_lowchol The lower Cholesky factor of the covariance matrix
#' @param x0 A vector (or matrix) of covariates at the prediction location(s)
#' @param c0 A vector (or matrix) of covariances between the prediction location(s) and the observed data
#'
#' @return The variance adjustment term that accounts for estimation of the
#'   latent random effects when predicting from a GLM-type model
#'
#' @noRd
get_wts_varw <- function(family, Xmat, y, w, size, dispersion, cov_lowchol, x0, c0) {
  # cov_lowchol is always a plain base matrix by the time this is called (see
  # get_pred_spglm), so base:: linear algebra avoids Matrix's S4 dispatch/
  # validity-check overhead -- this runs once per prediction row for the
  # local(distance/covariance) big-data approximation
  SigInv <- base::chol2inv(base::t(cov_lowchol)) # works on upchol
  SigInv_X <- SigInv %*% Xmat
  cov_betahat <- base::chol2inv(base::chol(base::crossprod(Xmat, SigInv_X)))
  wts_beta <- base::tcrossprod(cov_betahat, SigInv_X)
  # projection matrix that removes the fixed-effect (X) component of
  # SigInv, leaving only the part attributable to the latent random effects
  Ptheta <- SigInv - SigInv_X %*% wts_beta

  # d is the first derivative (score) of the log-likelihood w.r.t. the latent
  # predictor w; only needed conceptually for the gradient (commented below)
  d <- get_d(family, w, y, size, dispersion)
  # and then the gradient vector
  # g <-  d - Ptheta %*% w
  # Next, compute H
  # D is the second-derivative (curvature) contribution from the GLM likelihood
  D <- get_D(family, w, y, size, dispersion)
  # H is the Hessian of the penalized (Laplace) log-likelihood w.r.t. w;
  # -H (below) approximates the posterior precision of the latent effects
  H <- D - Ptheta
  mHInv <- base::solve(-H) # chol2inv(chol(Matrix::forceSymmetric(-H))) # solve(-H)

  if (is.vector(x0)) { # for length-one predicts result x0 c0 are vectors (how splm pred operates)
    # weight vector combining the fixed-effect contribution (x0) and the
    # kriging-type contribution (c0) to the linear predictor at the new location
    wts_pred <- x0 %*% wts_beta + c0 %*% SigInv - (c0 %*% SigInv_X) %*% wts_beta
    # quadratic form wts_pred' * mHInv * wts_pred gives the extra prediction
    # variance from having estimated (not known) the latent random effects
    var_adj <- as.numeric(wts_pred %*% base::tcrossprod(mHInv, wts_pred))
  } else { # this is to handle the matrix arguments for non-local predict calls with spglm
    if (NROW(x0) == 1) {
      wts_pred <- x0 %*% wts_beta + c0 %*% SigInv - (c0 %*% SigInv_X) %*% wts_beta
      var_adj <- as.numeric(wts_pred %*% base::tcrossprod(mHInv, wts_pred))
    } else {
      var_adj <- vapply(seq_len(NROW(x0)), function(x) { # this is so that only the diagonal of these products is returned
        x0_new <- x0[x, , drop = FALSE]
        c0_new <- c0[x, , drop = FALSE]
        wts_pred <- x0_new %*% wts_beta + c0_new %*% SigInv - (c0_new %*% SigInv_X) %*% wts_beta
        as.numeric(wts_pred %*% base::tcrossprod(mHInv, wts_pred))
      }, numeric(1))
    }
  }
  var_adj
}
