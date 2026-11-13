# Separate covariance contexts keep observed latent sampling wrt total covariance.
conditional_spatial_object <- function(object) {
  object$random <- NULL
  object$coefficients$randcov <- NULL
  object
}

get_conditional_random_design <- function(object, newdata) {
  observed <- if (inherits(object, c("spautor", "spgautor"))) object$data[object$observed_index, , drop = FALSE] else object$obdata
  get_randcov_design(observed, coef(object, type = "randcov"), object$partition_factor, newdata)
}


prepare_conditional_random <- function(object, design, covariance = NULL, spatial_factor = NULL) {
  Z <- design$observed
  n <- NROW(Z)
  active <- which(design$variance > 0 & design$info$conditioned)
  if (any(!is.finite(design$variance) | design$variance < 0)) stop("Invalid random-effect variances.", call. = FALSE)
  if (length(active)^2 > 1e8) warning(sprintf(
    "Conditional simulation requires factoring a covariance matrix for %s random-effect coefficients. This may result in exceedingly long computational times.",
    format(length(active), big.mark = ",", scientific = FALSE)
  ), call. = FALSE, immediate. = TRUE)
  spatial <- conditional_spatial_object(object)
  blocks <- list(seq_len(n))
  precision <- diag(1 / design$variance[active], length(active))
  prepared <- vector("list", length(blocks))
  for (k in seq_along(blocks)) {
    rows <- blocks[[k]]
    if (is.null(covariance)) {
      sub <- spatial
      sub$obdata <- spatial$obdata[rows, , drop = FALSE]
      S <- as.matrix(covmatrix(sub))
    } else S <- covariance[rows, rows, drop = FALSE]
    R <- if (is.null(spatial_factor)) conditional_local_factor(S, "Random-effect spatial block")$U else spatial_factor
    Sinv <- chol2inv(R)
    cols <- which(as.numeric(Matrix::colSums(abs(Z[rows, active, drop = FALSE]))) > 0)
    Zp <- as.matrix(Z[rows, active[cols], drop = FALSE])
    weights <- Sinv %*% Zp
    precision[cols, cols] <- precision[cols, cols, drop = FALSE] + crossprod(Zp, weights)
    prepared[[k]] <- list(rows = rows, cols = cols, weights = weights,
      factor = if (length(blocks) == 1L) R else NULL)
  }
  # One shared effect covariance: H = (Gamma^-1 + sum Zp' Sp^-1 Zp)^-1.
  H <- if (length(active)) chol2inv(conditional_local_factor(precision, "Random-effect precision")$U) else matrix(numeric(), 0, 0)
  L <- if (length(active)) t(conditional_local_factor(H, "Random-effect conditional covariance")$U) else matrix(numeric(), 0, 0)
  list(design = design, active = active, blocks = prepared, covariance = H, factor = L, spatial = spatial)
}

draw_conditional_random <- function(context, residual) {
  q <- length(context$design$variance)
  samples <- NCOL(residual)
  z <- matrix(rnorm(q * samples), q, samples)
  a <- sqrt(context$design$variance) * z
  if (length(context$active)) {
    rhs <- matrix(0, length(context$active), samples)
    for (block in context$blocks) {
      rhs[block$cols, ] <- rhs[block$cols, , drop = FALSE] +
        crossprod(block$weights, residual[block$rows, , drop = FALSE])
    }
    a[context$active, ] <- context$covariance %*% rhs + context$factor %*% z[context$active, , drop = FALSE]
  }
  a
}

conditional_random_output <- function(a, design) {
  rownames(a) <- rownames(design$info)
  attr(a, "effect_info") <- design$info
  a
}

check_conditional_random_output <- function(value) {
  if (is.null(value)) stop("Requested random-effect draws were not generated.", call. = FALSE)
  value
}

prepare_conditional_random_prediction <- function(context, newdata) {
  if (!NROW(newdata)) return(list())
  lower <- t(context$blocks[[1L]]$factor)
  covariance <- get_conditional_prediction_covariance(newdata, context$spatial)
  factor <- factor_conditional_prediction_covariance(covariance, lower)
  list(factor = factor, lower = lower)
}

draw_conditional_random_prediction <- function(context, prediction, residual, a) {
  m <- NROW(context$design$newdata)
  samples <- NCOL(residual)
  if (!m) return(matrix(numeric(), 0, samples))
  spatial <- residual - as.matrix(context$design$observed %*% a)
  whitened <- forwardsolve(prediction$lower, spatial)
  value <- unname(crossprod(prediction$factor$SqrtSigInv_c0, whitened) +
    prediction$factor$chol_cond_cov %*% matrix(rnorm(m * samples), m, samples))
  value + as.matrix(context$design$newdata %*% a)
}

get_conditional_random_lm <- function(object, newdata, Xnew, y, samples, design, random = FALSE) {
  X <- model.matrix(object)
  context <- prepare_conditional_random(object, design)
  prediction <- prepare_conditional_random_prediction(context, newdata)
  betahat <- coef(object)
  factor <- t(chol(vcov(object)))
  beta <- matrix(NA_real_, length(betahat), samples, dimnames = list(names(betahat), NULL))
  value <- matrix(NA_real_, NROW(Xnew), samples)
  effects <- if (random) matrix(NA_real_, length(design$variance), samples) else NULL
  chunk_size <- samples
  for (start in seq.int(1L, samples, by = chunk_size)) {
    cols <- seq.int(start, min(samples, start + chunk_size - 1L))
    beta[, cols] <- betahat + factor %*% matrix(rnorm(length(betahat) * length(cols)), length(betahat))
    residual <- y - X %*% beta[, cols, drop = FALSE]
    a <- draw_conditional_random(context, residual)
    value[, cols] <- Xnew %*% beta[, cols, drop = FALSE] +
      draw_conditional_random_prediction(context, prediction, residual, a)
    if (random) effects[, cols] <- a
  }
  list(beta = beta, newdata = value, random = if (random) conditional_random_output(effects, design) else NULL)
}

get_conditional_prediction_covariance <- function(newdata, object) {
  # Preserve new partition labels in covmatrix's prediction-prediction path.
  for (name in all.vars(object$partition_factor)) {
    if (is.character(newdata[[name]])) newdata[[name]] <- factor(newdata[[name]])
  }
  list(cross = covmatrix(object, newdata, cov_type = "obs.pred"),
    new = covmatrix(object, newdata, cov_type = "pred.pred"))
}

factor_conditional_prediction_covariance <- function(covariance, cov_lowchol_observed) {
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_observed, covariance$cross)
  cond_cov <- covariance$new - crossprod(SqrtSigInv_c0)
  list(SqrtSigInv_c0 = SqrtSigInv_c0, chol_cond_cov = t(chol(cond_cov)))
}
