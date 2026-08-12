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
#' \code{conditional.splm()}. Observed-data residuals against simulated draws
#' of beta play the role that a single fixed \code{base_val} plays in
#' \code{get_conditional_new_from_base()}. Because betahat's estimation
#' uncertainty is already propagated by that upstream simulation (drawing a
#' new beta and adding its trend back in is, by the law of total variance,
#' equivalent to adding an analytic fixed-effect-uncertainty term to a
#' fixed-betahat conditional covariance -- see the comments in
#' \code{conditional.splm()}), the ordinary kriging conditional covariance is
#' used here as-is, with no separate analytic correction. Adding one on top
#' of the simulated beta draws would double-count the same uncertainty.
#'
#' The conditional mean is linear in the per-draw residual \code{y_base - X_base
#' \%*\% beta_b}, so rather than solving/crossprod-ing an \code{n_base x
#' samples} residual matrix against the base covariance (repeating that
#' \code{O(n_base^2 * samples)}/\code{O(n_base * n_new * samples)} work for
#' every block, since \code{y_base}/\code{X_base}/\code{cov_lowchol_base} are
#' identical across blocks), the code solves against \code{y_base}
#' (\code{n_base x 1}) and \code{X_base} (\code{n_base x p}) once up front and
#' then recombines with each block's own \code{new_betahat} only
#' after the (block-specific, but \code{samples}-independent) crossprod with
#' \code{SqrtSigInv_c0}, an \code{O(n_new * p * samples)} recombination
#' instead, with \code{p} (the number of fixed effects) typically far smaller
#' than \code{samples}.
#'
#' @param newdata_list A list with elements \code{x0} (the newdata design
#'   matrix for this block; unused here, kept for a consistent calling
#'   convention with \code{\link{get_conditional_new_from_base_adjust_glm}()})
#'   and \code{newdata} (the newdata rows for this block).
#' @param object A fitted \code{splm} model object.
#' @param SqrtSigInv_y \code{forwardsolve(cov_lowchol_base, y_base)}, computed
#'   once by the caller (shared across every block).
#' @param SqrtSigInv_X \code{forwardsolve(cov_lowchol_base, X_base)}, computed
#'   once by the caller (shared across every block).
#' @param new_betahat A matrix of simulated beta draws, one column per
#'   simulation.
#' @param cov_lowchol_base The lower triangular Cholesky factor of the base
#'   locations' covariance matrix.
#' @param samples The number of simulations.
#'
#' @return A matrix of simulated residuals at this block's \code{newdata}
#'   rows, one column per simulation.
#'
#' @noRd
get_conditional_new_from_base_adjust <- function(newdata_list, object, SqrtSigInv_y, SqrtSigInv_X, new_betahat, cov_lowchol_base, samples) {

  newdata <- newdata_list$newdata

  newdata_n <- NROW(newdata)

  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")


  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, cov_base_new)

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

  # conditional mean: Sigma_{new,base} Sigma_base^-1 (y_base - X_base %*%
  # beta_b) == crossprod(SqrtSigInv_c0, SqrtSigInv_y) - crossprod(SqrtSigInv_c0,
  # SqrtSigInv_X) %*% beta_b
  cond_mu_y <- crossprod(SqrtSigInv_c0, SqrtSigInv_y)
  cond_mu_X <- crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  cond_mu <- as.numeric(cond_mu_y) - cond_mu_X %*% new_betahat
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
#' As in \code{\link{get_conditional_new_from_base_adjust}()}, the conditional
#' mean is linear in the per-draw residual \code{w_base - X_base \%*\% beta_b},
#' so the code solves \code{cov_lowchol_base} against \code{w_base} and
#' \code{X_base} once (the latter, \code{SqrtSigInv_X}, is already computed by
#' \code{conditional.spglm()} for \code{var_adj} and simply reused here) and
#' this function recombines with each block's own \code{new_betahat} only
#' after crossprod-ing with \code{SqrtSigInv_c0}.
#'
#' @param newdata_list A list with elements \code{x0} (the newdata design
#'   matrix for this block) and \code{newdata} (the newdata rows for this
#'   block).
#' @param object A fitted \code{spglm} model object.
#' @param SqrtSigInv_w \code{forwardsolve(cov_lowchol_base, w_base)}, computed
#'   once by the caller (shared across every block).
#' @param SqrtSigInv_X \code{forwardsolve(cov_lowchol_base, X_base)}, computed
#'   once by the caller for \code{var_adj} and reused here (shared across
#'   every block).
#' @param new_betahat A matrix of simulated beta draws, one column per
#'   simulation.
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
get_conditional_new_from_base_adjust_glm <- function(newdata_list, object, SqrtSigInv_w, SqrtSigInv_X, new_betahat, cov_lowchol_base, samples, SigInv, SigInv_X, wts_beta, cov_lowchol_mH) {

  x0 <- newdata_list$x0
  newdata <- newdata_list$newdata

  newdata_n <- NROW(newdata)
  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred")
  cov_new <- covmatrix(object, newdata, cov_type = "pred.pred")


  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, cov_base_new)

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

  # conditional mean, split into a fixed (not-per-sample) piece against
  # w_base and a cheap p-column piece against X_base
  cond_mu_w <- crossprod(SqrtSigInv_c0, SqrtSigInv_w)
  cond_mu_X <- crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  cond_mu <- as.numeric(cond_mu_w) - cond_mu_X %*% new_betahat
  new_val <- new_val + cond_mu
}
