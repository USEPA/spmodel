# Shared LM/GLM block draws; residuals have already been whitened on the base.
draw_conditional_lowrank_blocks <- function(object, newdata, Xnew, beta, residual, cov_lowchol_base, local_list, effects = NULL, Znew = NULL) {
  samples <- NCOL(beta)
  chunk_size <- get_conditional_chunk_size(local_list)
  new_val <- matrix(NA_real_, NROW(newdata), samples, dimnames = list(rownames(Xnew), NULL))
  index <- if (!NROW(newdata)) list() else if (local_list$method_new == "all") list(seq_len(NROW(newdata))) else local_list$index$new
  if (isTRUE(local_list$parallel)) {
    cl <- parallel::makeCluster(local_list$ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  width <- if (local_list$parallel) local_list$ncores else 1L
  for (first in if (length(index)) seq.int(1L, length(index), by = width) else integer()) {
    batch <- seq.int(first, min(length(index), first + width - 1L))
    if (local_list$parallel) {
      # Bound worker payloads and avoid serializing the fitted formula environment.
      covariance <- lapply(index[batch], function(rows) {
        get_conditional_glm_block_cov(newdata[rows, , drop = FALSE], object)
      })
      blocks <- parallel::parLapply(cl, covariance, factor_conditional_glm_block,
        cov_lowchol_base)
      rm(covariance)
    } else {
      blocks <- list(get_conditional_glm_block(newdata[index[[first]], , drop = FALSE],
        object, cov_lowchol_base))
    }
    for (j in seq_along(batch)) {
      rows <- index[[batch[j]]]
      block <- blocks[[j]]
      for (start in seq.int(1L, samples, by = chunk_size)) {
        columns <- seq.int(start, min(samples, start + chunk_size - 1L))
        new_val[rows, columns] <- Xnew[rows, , drop = FALSE] %*% beta[, columns, drop = FALSE] +
          crossprod(block$SqrtSigInv_c0, residual[, columns, drop = FALSE]) + block$chol_cond_cov %*%
          matrix(rnorm(length(rows) * length(columns)), length(rows), length(columns))
        if (!is.null(effects)) new_val[rows, columns] <- new_val[rows, columns, drop = FALSE] +
          as.matrix(Znew[rows, , drop = FALSE] %*% effects[, columns, drop = FALSE])
      }
    }
    rm(blocks, block)
  }
  new_val
}
