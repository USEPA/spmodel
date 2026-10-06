# Prepare the complete observed latent sequence. Response neighbors K may be
# undrawn; latent neighbors N must precede i. Two separate steps
get_conditional_vecchia_glm_context <- function(object, newdata, Xnew, local_list) {
  X <- model.matrix(object)
  n <- NROW(X)
  eta <- fitted(object, type = "link")
  w <- as.numeric(w_offset_free(eta, model.offset(model.frame(object))))
  D <- diag(get_D(object$family, eta, object$y, object$size, as.vector(coef(object, type = "dispersion"))))
  if (any(!is.finite(w)) || any(!is.finite(D))) stop("Observed latent values and response curvature must be finite.", call. = FALSE)
  covariance <- get_conditional_vecchia_covariance(object, object$obdata)
  ord <- conditional_vecchia_order(covariance$coords, local_list$ordering)
  operators <- vector("list", n)
  all_sites <- seq_len(n)
  for (k in seq_len(n)) {
    i <- ord[k]
    # Compute scores once over the observed sites. N is restricted to predecessors;
    # K may include undrawn sites. The current response uses one slot in K's budget,
    # and overlap between N and K is not refilled. Thus |A| <= 2 * size.
    scores <- if (local_list$method != "all") conditional_vecchia_scores(covariance, i, local_list$method) else NULL
    N <- conditional_vecchia_neighbors(covariance, i, ord[seq_len(k - 1L)], local_list$size, local_list$method, scores)
    K <- c(i, conditional_vecchia_neighbors(covariance, i, all_sites, local_list$size - 1L, local_list$method, scores))
    A <- unique(c(i, N, K))
    label <- paste("Observed site", i, "neighborhood")
    # Step 1 constructs a proper Gaussian, centered at fitted w_A plus M_A delta.
    gaussian <- prepare_conditional_local_gaussian(conditional_vecchia_covariance(covariance, A),
      X[A, , drop = FALSE], w[A], D[A], label)
    # Step 2 conditions that same Gaussian on N and caches only the target.
    op <- prepare_conditional_latent_operator(gaussian, 1L, match(N, A), label)
    op$i <- i
    op$neighbors <- N
    operators[[k]] <- op
  }
  prediction <- if (NROW(newdata)) prepare_conditional_vecchia(object, newdata, local_list) else NULL
  list(operators = operators, order = ord, prediction = prediction, X = X,
    betahat = coef(object), coefficient_factor = t(chol(vcov(object))))
}

# Draw shared coefficients, observed latent values, and prediction residuals.
# Optional latent output retains the field used for prediction.
get_conditional_vecchia_glm <- function(object, newdata, Xnew, local_list, samples, latent = FALSE) {
  chunk_size <- get_conditional_glm_chunk_size(local_list)
  context <- get_conditional_vecchia_glm_context(object, newdata, Xnew, local_list)
  n <- NROW(context$X)
  m <- NROW(Xnew)
  p <- length(context$betahat)
  beta <- matrix(NA_real_, p, samples, dimnames = list(names(context$betahat), NULL))
  new_val <- matrix(NA_real_, m, samples, dimnames = list(rownames(Xnew), NULL))
  observed <- if (latent) matrix(NA_real_, n, samples, dimnames = list(rownames(context$X), NULL)) else NULL
  # Matrix preparation is outside this loop. Chunking bounds working storage;
  # every column gets one shared coefficient draw for all observed/prediction sites.
  for (start in seq.int(1L, samples, by = chunk_size)) {
    columns <- seq.int(start, min(samples, start + chunk_size - 1L))
    count <- length(columns)
    delta <- context$coefficient_factor %*% matrix(rnorm(p * count), p, count)
    beta[, columns] <- sweep(delta, 1L, context$betahat, "+")
    w <- matrix(NA_real_, n, count)
    # Apply w_i = c_i + t_i delta + g_i w_N + sqrt(v_i) e_i to all columns.
    # Only predecessors in N have been drawn; the remaining sites of A contributed
    # to the prepared Gaussian but are marginalized, not separately sampled here.
    for (op in context$operators) {
      mu <- op$intercept + as.numeric(crossprod(op$coefficient, delta))
      if (length(op$neighbors)) mu <- mu + as.numeric(crossprod(op$weights, w[op$neighbors, , drop = FALSE]))
      w[op$i, ] <- mu + sqrt(op$variance) * rnorm(count)
    }
    if (latent) observed[, columns] <- w
    if (m) {
      # Once w is drawn, w - X beta is known within each simulation. Reuse the
      # spatial residual conditionals, then restore the same coefficient trend.
      residual <- draw_conditional_vecchia(context$prediction, w - context$X %*% beta[, columns, drop = FALSE], count)
      new_val[, columns] <- Xnew %*% beta[, columns, drop = FALSE] + residual
    }
  }
  list(beta = beta, newdata = new_val, latent = observed)
}
