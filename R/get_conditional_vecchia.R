# Prepare spatial conditionals once. Every observed site and earlier prediction
# is a candidate; only the selected neighborhood covariance is factored.
# The weights and variances depend on the fitted covariance and locations, not on
# simulated residual values, so all samples and GLM chunks share these operators.
# Covariance ranking includes random effects and partition eligibility; finite
# neighborhoods still approximate the information carried by global random effects.
prepare_conditional_vecchia <- function(object, newdata, local_list) {
  n <- NROW(object$obdata)
  m <- NROW(newdata)
  columns <- unique(c(object$xcoord, object$ycoord, all.vars(object$random), all.vars(object$partition_factor)))
  data <- rbind(object$obdata[, columns, drop = FALSE], newdata[local_list$order, columns, drop = FALSE])
  context <- get_conditional_vecchia_covariance(object, data)
  operators <- vector("list", m)
  for (k in seq_len(m)) {
    i <- n + k
    neighbors <- conditional_vecchia_neighbors(context, i, seq_len(i - 1L), local_list$size, local_list$method)
    covariance <- conditional_vecchia_covariance(context, c(i, neighbors))
    op <- prepare_conditional_local_weights(covariance, 1L, seq_along(neighbors) + 1L,
      paste("Prediction site", local_list$order[k]))
    op$neighbors <- neighbors
    operators[[k]] <- op
  }
  list(operators = operators, order = local_list$order, n = n)
}

# Apply cached operators across a simulation chunk. Each realized residual is
# fixed when conditioning later sites; no extra latent uncertainty is added.
draw_conditional_vecchia <- function(context, base_val, samples) {
  m <- length(context$order)
  pool <- matrix(NA_real_, context$n + m, samples)
  pool[seq_len(context$n), ] <- base_val
  Z <- matrix(rnorm(m * samples), m, samples)
  # Each row uses already realized residuals with g_i = Sigma_iN Sigma_NN^-1
  # and v_i = Sigma_ii - g_i Sigma_Ni. The new innovation is independent of the
  # preceding ones, and the saved result becomes available to subsequent sites.
  for (k in seq_len(m)) {
    op <- context$operators[[k]]
    pool[context$n + k, ] <- as.numeric(crossprod(op$weights, pool[op$neighbors, , drop = FALSE])) +
      sqrt(op$variance) * Z[k, ]
  }
  value <- matrix(NA_real_, m, samples)
  value[context$order, ] <- pool[context$n + seq_len(m), , drop = FALSE]
  value
}

# Linear models prepare once per call; GLMs reuse the preparation across chunks.
get_conditional_vecchia <- function(object, newdata, base_val, local_list, samples) {
  draw_conditional_vecchia(prepare_conditional_vecchia(object, newdata, local_list), base_val, samples)
}
