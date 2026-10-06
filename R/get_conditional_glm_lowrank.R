#' Prepare the base-sized GLM conditional approximation
#' @param object A fitted spglm object.
#' @param index Observed rows selected for the simulation base.
#' @return Whitened conditional factors and a base-sized covariance context.
#' @noRd
get_conditional_glm_base <- function(object, index) {
  X <- model.matrix(object)[index, , drop = FALSE]
  eta <- fitted(object, type = "link")[index]
  w <- w_offset_free(fitted(object, type = "link"), model.offset(model.frame(object)))[index]
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
    wts_latent = wts_latent, betahat = coef(object))
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
get_conditional_glm_lowrank <- function(object, newdata, Xnew, local_list, samples, latent = FALSE) {
  chunk_size <- get_conditional_glm_chunk_size(local_list)
  base <- get_conditional_glm_base(object, local_list$index$base)
  index <- if (!NROW(newdata)) list() else if (local_list$method_new == "all") list(seq_len(NROW(newdata))) else local_list$index$new
  if (local_list$parallel) {
    cl <- parallel::makeCluster(local_list$ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  n_base <- NROW(base$cov_lowchol_base)
  p <- length(base$betahat)
  beta <- matrix(NA_real_, p, samples, dimnames = list(names(base$betahat), NULL))
  new_val <- matrix(NA_real_, NROW(newdata), samples, dimnames = list(rownames(Xnew), NULL))
  residual <- matrix(NA_real_, n_base, samples)
  residual_mean <- base$SqrtSigInv_w - as.numeric(base$SqrtSigInv_X %*% base$betahat)
  residual_weights <- base$wts_latent - base$SqrtSigInv_X
  for (start in seq.int(1L, samples, by = chunk_size)) {
    columns <- seq.int(start, min(samples, start + chunk_size - 1L))
    n_draw <- length(columns)
    delta <- base$cov_betahat_lowchol %*% matrix(rnorm(p * n_draw), p, n_draw)
    beta[, columns] <- sweep(delta, 1L, base$betahat, "+")
    z <- backsolve(t(base$cov_lowchol_mH), matrix(rnorm(n_base * n_draw), n_base, n_draw))
    residual[, columns] <- sweep(residual_weights %*% delta + z, 1L, residual_mean, "+")
  }
  width <- if (local_list$parallel) local_list$ncores else 1L
  for (first in if (length(index)) seq.int(1L, length(index), by = width) else integer()) {
    batch <- seq.int(first, min(length(index), first + width - 1L))
    if (local_list$parallel) {
      # Bound worker payloads and avoid serializing the fitted formula environment.
      covariance <- lapply(index[batch], function(rows) {
        get_conditional_glm_block_cov(newdata[rows, , drop = FALSE], base$covariance)
      })
      blocks <- parallel::parLapply(cl, covariance, factor_conditional_glm_block,
        base$cov_lowchol_base)
      rm(covariance)
    } else {
      blocks <- list(get_conditional_glm_block(newdata[index[[first]], , drop = FALSE],
        base$covariance, base$cov_lowchol_base))
    }
    for (j in seq_along(batch)) {
      rows <- index[[batch[j]]]
      block <- blocks[[j]]
      for (start in seq.int(1L, samples, by = chunk_size)) {
        columns <- seq.int(start, min(samples, start + chunk_size - 1L))
        new_val[rows, columns] <- Xnew[rows, , drop = FALSE] %*% beta[, columns, drop = FALSE] +
          crossprod(block$SqrtSigInv_c0, residual[, columns, drop = FALSE]) + block$chol_cond_cov %*%
          matrix(rnorm(length(rows) * length(columns)), length(rows), length(columns))
      }
    }
    rm(blocks, block)
  }
  observed <- NULL
  if (latent) {
    # Recover the base draw from its whitened residual.
    rows <- local_list$index$base
    observed <- matrix(NA_real_, NROW(model.matrix(object)), samples)
    observed[rows, ] <- base$cov_lowchol_base %*% residual + model.matrix(object)[rows, , drop = FALSE] %*% beta
  }
  list(beta = beta, newdata = new_val, latent = observed)
}

#' Resolve the draw chunk size for local GLM simulation
#' @param local_list Resolved local simulation settings.
#' @return A positive integer chunk size.
#' @noRd
get_conditional_glm_chunk_size <- function(local_list) {
  chunk_size <- local_list$chunk_size
  if (is.null(chunk_size)) chunk_size <- 1000L
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || is.na(chunk_size) ||
      !is.finite(chunk_size) || chunk_size < 1 || chunk_size != floor(chunk_size)) {
    stop("local$chunk_size must be a positive integer.", call. = FALSE)
  }
  chunk_size
}
