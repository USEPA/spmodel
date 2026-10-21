#' Prepare the base-sized GLM conditional approximation
#' @param object A fitted spglm object.
#' @param index Observed rows selected for the simulation base.
#' @return Whitened conditional factors and a base-sized covariance context.
#' @noRd
get_conditional_glm_base <- function(object, index, center = get_latent_center(object)) {
  X <- model.matrix(object)[index, , drop = FALSE]
  eta <- fitted(object, type = "link")[index]
  w <- w_offset_free(center$eta, model.offset(model.frame(object)))[index]
  covariance <- object[c("obdata", "coefficients", "random", "partition_factor",
    "anisotropy", "xcoord", "ycoord", "dim_coords", "diagtol")]
  class(covariance) <- class(object)
  covariance$obdata <- object$obdata[index, , drop = FALSE]
  cov_lowchol_base <- t(chol(covmatrix(covariance)))
  SqrtSigInv_X <- forwardsolve(cov_lowchol_base, X)
  D <- get_D(object$family, eta, object$y[index], object$size[index],
    as.vector(coef(object, type = "dispersion")))
  # Congruence by the base factor avoids forming the spatial precision.
  mH <- diag(length(index)) - crossprod(cov_lowchol_base, D %*% cov_lowchol_base)
  cov_lowchol_mH <- tryCatch(t(chol(as.matrix(mH))), error = function(e) {
    stop("The simulation base latent precision is not positive definite; check the fitted model or select a different base.", call. = FALSE)
  })
  cov_betahat <- vcov(object)
  wts_latent <- backsolve(t(cov_lowchol_mH), forwardsolve(cov_lowchol_mH, SqrtSigInv_X))
  list(covariance = covariance, cov_lowchol_base = cov_lowchol_base,
    cov_lowchol_mH = cov_lowchol_mH, cov_betahat_lowchol = t(chol(cov_betahat)),
    SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_w = forwardsolve(cov_lowchol_base, w),
    wts_latent = wts_latent, betahat = center$betahat)
}

#' Prepare one low-rank GLM prediction block
#' @param newdata Prediction rows for this block.
#' @param object Base-sized covariance context.
#' @param cov_lowchol_base Lower base covariance factor.
#' @return Cached cross-covariance operator and conditional error factor.
#' @noRd
get_conditional_glm_block <- function(newdata, object, cov_lowchol_base) {
  covariance <- get_conditional_glm_block_cov(newdata, object)
  factor_conditional_glm_block(covariance, cov_lowchol_base)
}

#' @rdname get_conditional_glm_block
#' @noRd
get_conditional_glm_block_cov <- function(newdata, object) {
  # Preserve new partition labels in covmatrix's prediction-prediction path.
  for (name in all.vars(object$partition_factor)) {
    if (is.character(newdata[[name]])) newdata[[name]] <- factor(newdata[[name]])
  }
  list(cross = covmatrix(object, newdata, cov_type = "obs.pred"),
    new = covmatrix(object, newdata, cov_type = "pred.pred"))
}

#' @rdname get_conditional_glm_block
#' @param covariance Numeric cross-covariance and prediction covariance matrices.
#' @noRd
factor_conditional_glm_block <- function(covariance, cov_lowchol_base) {
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, covariance$cross)
  cond_cov <- covariance$new - crossprod(SqrtSigInv_c0)
  list(SqrtSigInv_c0 = SqrtSigInv_c0, chol_cond_cov = t(chol(cond_cov)))
}

#' Draw low-rank GLM simulations with cached factors
#' @param object A fitted spglm object.
#' @param newdata Processed prediction data.
#' @param Xnew Aligned prediction design matrix.
#' @param local_list Resolved low-rank settings.
#' @param samples Number of draws.
#' @return Matrices of coefficients and offset-free prediction latent values.
#' @noRd
get_conditional_glm_lowrank <- function(object, newdata, Xnew, local_list, samples, latent = FALSE, random_design = NULL, center = get_latent_center(object)) {
  if (!is.null(object$random) && is.null(random_design)) random_design <- get_conditional_random_design(object, newdata, local_list$index$base)
  chunk_size <- get_conditional_chunk_size(local_list)
  base <- get_conditional_glm_base(object, local_list$index$base, center)
  effect <- if (!is.null(object$random)) {
    sub <- object
    sub$obdata <- object$obdata[local_list$index$base, , drop = FALSE]
    prepare_conditional_random(sub, random_design)
  } else NULL
  n_base <- NROW(base$cov_lowchol_base)
  p <- length(base$betahat)
  beta <- matrix(NA_real_, p, samples, dimnames = list(names(base$betahat), NULL))
  residual <- matrix(NA_real_, n_base, samples)
  residual_mean <- base$SqrtSigInv_w - as.numeric(base$SqrtSigInv_X %*% base$betahat)
  residual_weights <- base$wts_latent - base$SqrtSigInv_X
  observed <- if (latent) matrix(NA_real_, NROW(model.matrix(object)), samples) else NULL
  effects <- if (!is.null(effect)) matrix(NA_real_, length(random_design$variance), samples) else NULL
  for (start in seq.int(1L, samples, by = chunk_size)) {
    columns <- seq.int(start, min(samples, start + chunk_size - 1L))
    n_draw <- length(columns)
    delta <- base$cov_betahat_lowchol %*% matrix(rnorm(p * n_draw), p, n_draw)
    beta[, columns] <- sweep(delta, 1L, base$betahat, "+")
    z <- backsolve(t(base$cov_lowchol_mH), matrix(rnorm(n_base * n_draw), n_base, n_draw))
    residual[, columns] <- sweep(residual_weights %*% delta + z, 1L, residual_mean, "+")
    if (latent) observed[local_list$index$base, columns] <- base$cov_lowchol_base %*% residual[, columns, drop = FALSE] +
      model.matrix(object)[local_list$index$base, , drop = FALSE] %*% beta[, columns, drop = FALSE]
    if (!is.null(effect)) {
      natural <- base$cov_lowchol_base %*% residual[, columns, drop = FALSE]
      a <- draw_conditional_random(effect, natural)
      effects[, columns] <- a
      residual[, columns] <- forwardsolve(t(effect$blocks[[1L]]$factor), natural - as.matrix(random_design$observed %*% a))
    }
  }
  if (!is.null(effect)) {
    base$covariance <- effect$spatial
    base$cov_lowchol_base <- t(effect$blocks[[1L]]$factor)
  }
  new_val <- draw_conditional_lowrank_blocks(base$covariance, newdata, Xnew, beta, residual,
    base$cov_lowchol_base, local_list, effects, if (!is.null(effect)) random_design$newdata else NULL)
  list(beta = beta, newdata = new_val, latent = observed)
}
