#' Simulate one block of a base-and-block conditional simulation
#'
#' Given a Gaussian process already simulated at a "base" set of locations
#' (\code{base_val}), draws values at a new set of locations (\code{newdata})
#' from their conditional distribution given the base draws, treating the
#' mean (fixed effect trend) as known/already removed from \code{base_val}.
#' This is the workhorse behind \code{\link{sprnorm}()}'s big data
#' approximation: simulating the full field jointly is
#' \eqn{O((n_{base} + n_{new})^3)}, while simulating the base sample once and
#' then each block conditionally is much cheaper, at the cost of treating
#' distinct blocks as conditionally independent given the base sample. See
#' \code{\link{get_conditional_new_from_base_adjust}()} for the variant used
#' by \code{\link{conditional}()}, which additionally propagates fixed effect
#' (beta) uncertainty into the simulated values.
#'
#' @param newdata The new locations' data (only used for its row count).
#' @param object A fitted model object.
#' @param base_val A matrix of simulated (mean-zero) values at the base
#'   locations, one column per simulation.
#' @param cov_lowchol_base The lower triangular Cholesky factor of the base
#'   locations' covariance matrix.
#' @param samples The number of simulations (i.e. columns of \code{base_val}).
#'
#' @return A matrix of simulated values at \code{newdata}, one column per
#'   simulation.
#'
#' @noRd
get_conditional_new_from_base <- function(newdata, object, base_val, cov_lowchol_base, samples) {

  newdata_n <- NROW(newdata)
  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")

  # standard Gaussian conditioning (kriging) formulas, computed via the
  # Cholesky factor rather than an explicit inverse of the base covariance:
  # SqrtSigInv_c0 = Sigma_base^(-1/2) %*% Sigma_{base,new}
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, cov_base_new)
  SqrtSigInv_base_val <- forwardsolve(cov_lowchol_base, base_val)

  # conditional covariance: Sigma_new - Sigma_{new,base} Sigma_base^-1 Sigma_{base,new}
  cond_cov <- cov_new - crossprod(SqrtSigInv_c0, SqrtSigInv_c0)
  chol_cond_cov <- t(chol(cond_cov))
  new_val <- vapply(seq_len(samples), function(x) as.numeric(chol_cond_cov %*% rnorm(newdata_n)), numeric(newdata_n))

  # conditional mean: Sigma_{new,base} Sigma_base^-1 base_val, added onto the
  # mean-zero draw from cond_cov above
  cond_mu <- crossprod(SqrtSigInv_c0, SqrtSigInv_base_val)
  new_val <- new_val + cond_mu
}

#' Simulate one block of a base-and-block conditional simulation, for
#' \code{splm} models
#'
#' Variant of \code{\link{get_conditional_new_from_base}()} used by
#' \code{conditional.splm()}. \code{base_val} there holds observed-data
#' residuals computed against simulated draws of beta rather than a single
#' fixed betahat. Because betahat's estimation uncertainty is already
#' propagated by that upstream simulation (drawing a new beta and adding its
#' trend back in is, by the law of total variance, equivalent to adding an
#' analytic fixed-effect-uncertainty term to a fixed-betahat conditional
#' covariance -- see the comments in \code{conditional.splm()}), the ordinary
#' kriging conditional covariance is used here as-is, with no separate
#' analytic correction. Adding one on top of the simulated beta draws would
#' double-count the same uncertainty.
#'
#' @param newdata_list A list with elements \code{x0} (the newdata design
#'   matrix for this block; unused here, kept for a consistent calling
#'   convention with \code{\link{get_conditional_new_from_base_adjust_glm}()})
#'   and \code{newdata} (the newdata rows for this block).
#' @param object A fitted \code{splm} model object.
#' @param base_val A matrix of simulated base-sample residuals (one column
#'   per simulated beta draw; see Details above).
#' @param cov_lowchol_base The lower triangular Cholesky factor of the base
#'   locations' covariance matrix.
#' @param samples The number of simulations.
#'
#' @return A matrix of simulated residuals at this block's \code{newdata}
#'   rows, one column per simulation.
#'
#' @noRd
get_conditional_new_from_base_adjust <- function(newdata_list, object, base_val, cov_lowchol_base, samples) {

  newdata <- newdata_list$newdata

  newdata_n <- NROW(newdata)

  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")


  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, cov_base_new)
  SqrtSigInv_base_val <- forwardsolve(cov_lowchol_base, base_val)

  # ordinary kriging conditional covariance -- betahat uncertainty is already
  # supplied by the caller's simulated beta draws, so it is not added here
  cond_cov <- cov_new - crossprod(SqrtSigInv_c0, SqrtSigInv_c0)

  spcov_val <- coef(object, type = "spcov")
  # pure nugget: cond_cov is already diagonal (no spatial dependence to
  # condition on), so take element-wise sqrt() instead of a full chol()
  if (spcov_val[["de"]] == 0 && is.null(coef(object, type = "randcov"))) {
    chol_cond_cov <- Matrix::Diagonal(NROW(cond_cov))
    diag(chol_cond_cov) <- sqrt(diag(chol_cond_cov))
  } else {
    chol_cond_cov <- t(chol(cond_cov))
  }

  new_val <- vapply(seq_len(samples), function(x) as.numeric(chol_cond_cov %*% rnorm(newdata_n)), numeric(newdata_n))

  cond_mu <- crossprod(SqrtSigInv_c0, SqrtSigInv_base_val)
  new_val <- new_val + cond_mu
}

#' Simulate one block of a base-and-block conditional simulation, adjusted
#' for latent-process (\code{w}) estimation uncertainty (\code{spglm} models)
#'
#' GLM analog of \code{\link{get_conditional_new_from_base_adjust}()}, used
#' by \code{conditional.spglm()}. As in the Gaussian case, fixed-effect
#' (beta) uncertainty is supplied purely by the caller's simulated beta
#' draws, so no analytic fixed-effect correction is added here. The
#' link-scale latent process \code{w}, however, is held fixed at its fitted
#' value by the caller (not simulated) -- its own Laplace-approximate
#' posterior uncertainty (\code{cov_lowchol_mH}) is instead supplied
#' analytically here via \code{var_adj}, computed from the prediction weights
#' implied by that approximation. Simulating a new \code{w} upstream and
#' including \code{var_adj} here would double-count the same uncertainty.
#'
#' @param newdata_list A list with elements \code{x0} (the newdata design
#'   matrix for this block) and \code{newdata} (the newdata rows for this
#'   block).
#' @param object A fitted \code{spglm} model object.
#' @param base_val A matrix of simulated base-sample link-scale residuals
#'   (one column per simulated beta draw).
#' @param cov_lowchol_base The lower triangular Cholesky factor of the base
#'   locations' covariance matrix.
#' @param samples The number of simulations.
#' @param SigInv The precision matrix of the base locations' covariance
#'   matrix.
#' @param SigInv_X \code{Sigma_base^-1 \%*\% X}.
#' @param wts_beta \code{cov_betahat \%*\% t(SigInv_X)}, prediction weights
#'   for the fixed effect contribution.
#' @param cov_lowchol_mH The lower triangular Cholesky factor of the negative
#'   Hessian of the joint log-likelihood for \code{w} (see
#'   \code{conditional.spglm()}), used to weight the \code{var_adj}
#'   adjustment below.
#'
#' @return A matrix of simulated link-scale residuals at this block's
#'   \code{newdata} rows, one column per simulation.
#'
#' @noRd
get_conditional_new_from_base_adjust_glm <- function(newdata_list, object, base_val, cov_lowchol_base, samples, SigInv, SigInv_X, wts_beta, cov_lowchol_mH) {

  x0 <- newdata_list$x0
  newdata <- newdata_list$newdata

  newdata_n <- NROW(newdata)
  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")


  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, cov_base_new)
  SqrtSigInv_base_val <- forwardsolve(cov_lowchol_base, base_val)

  # ordinary kriging conditional covariance -- betahat uncertainty is already
  # supplied by the caller's simulated beta draws, so it is not added here
  cond_cov <- cov_new - crossprod(SqrtSigInv_c0, SqrtSigInv_c0)

  # prediction weights for the new locations under the Laplace
  # approximation's linear predictor for w (same structure as the universal
  # kriging weights used elsewhere in the package, combining a fixed effect
  # term and a covariance-based term)
  c0 <- t(cov_base_new)
  wts_pred <- x0 %*% wts_beta + c0 %*% SigInv - (c0 %*% SigInv_X) %*% wts_beta
  wts_pred <- t(wts_pred)
  # project those weights through the Laplace posterior precision's Cholesky
  # factor to get the additional predictive variance contributed by not
  # knowing w exactly (only its Laplace-approximate posterior)
  SqrtmHInv_wts_pred <- forwardsolve(cov_lowchol_mH, wts_pred)
  var_adj <- crossprod(SqrtmHInv_wts_pred, SqrtmHInv_wts_pred)

  # fold the linearization-uncertainty adjustment into the conditional
  # covariance before factoring it, so the extra uncertainty from only
  # knowing w up to its Laplace-approximate posterior actually widens the
  # simulated draws below
  cond_cov <- var_adj + cond_cov

  chol_cond_cov <- t(chol(cond_cov))
  new_val <- vapply(seq_len(samples), function(x) as.numeric(chol_cond_cov %*% rnorm(newdata_n)), numeric(newdata_n))

  cond_mu <- crossprod(SqrtSigInv_c0, SqrtSigInv_base_val)
  new_val <- new_val + cond_mu
}
