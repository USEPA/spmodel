# factor conditional matrix
conditional_local_factor <- function(A, label) {
  A <- as.matrix(A)
  n <- NROW(A)
  fail <- function(reason) stop(paste0(label, ": ", reason), call. = FALSE)
  if (!n || NCOL(A) != n || any(!is.finite(A))) fail("matrix must be finite and square.")
  if (any(diag(A) <= 0)) fail("matrix is not positive definite.")
  U <- tryCatch(base::chol(A), error = function(e) fail("matrix is not positive definite."))
  list(U = U)
}

# Solve using triangular factor
conditional_local_solve <- function(factor, rhs) {
  rhs <- as.matrix(rhs)
  backsolve(factor$U, forwardsolve(t(factor$U), rhs))
}

# Step 1: form the fitted-centered local Gaussian. D uses fitted link values;
# spatial w excludes offsets.
prepare_conditional_local_gaussian <- function(Sigma, X, w, D, label = "Local latent Gaussian") {
  spatial <- conditional_local_factor(Sigma, paste(label, "spatial covariance"))
  # chol2inv reuses the factor and avoids solving against an identity matrix.
  SigInv <- base::chol2inv(spatial$U)
  precision <- SigInv
  diag(precision) <- diag(precision) - D
  factor <- conditional_local_factor(precision, paste(label, "latent precision"))
  V <- base::chol2inv(factor$U)
  # Retain the precision factor
  list(mean = as.numeric(w), M = V %*% (SigInv %*% X), V = V, factor = factor)
}

# Condition on neighbor rows
prepare_conditional_local_weights <- function(V, i, N, label, precision_factor = NULL) {
  if (!length(N)) return(list(weights = numeric(), variance = V[i, i]))
  factor <- conditional_local_factor(V[N, N, drop = FALSE], paste(label, "neighbor covariance"))
  g <- as.numeric(conditional_local_solve(factor, V[N, i, drop = FALSE]))
  explained <- sum(g * V[N, i])
  variance <- V[i, i] - explained
  # Allow scaling of rounding error
  tol <- 100 * NROW(V) * .Machine$double.eps * max(abs(V[i, i]), abs(explained))
  if (!is.finite(variance) || variance < -tol) stop(paste0(label, ": conditional variance is negative or nonfinite."), call. = FALSE)
  if (variance <= sqrt(.Machine$double.eps) * max(abs(V[i, i]), abs(explained))) {
    if (is.null(precision_factor)) {
      full <- conditional_local_factor(V, paste(label, "full covariance"))
      K <- t(full$U)
    } else {
      K <- backsolve(precision_factor$U, diag(NROW(V)))
    }
    projection <- qr(t(K[N, , drop = FALSE]), tol = 100 * NROW(V) * .Machine$double.eps)
    if (projection$rank != length(N)) stop(paste0(label, ": numerical neighbor rank is deficient."), call. = FALSE)
    variance <- sum(qr.qty(projection, K[i, ])[-seq_along(N)]^2)
  }
  if (!is.finite(variance) || variance <= 0) stop(paste0(label, ": conditional variance cannot be resolved at numerical precision."), call. = FALSE)
  list(weights = g, variance = variance)
}

# Step 2: cache the conditional intercept and weights for the current site.
# This deliberately keeps the fitted center, with no gradient correction c_A.
prepare_conditional_latent_operator <- function(gaussian, i, N, label) {
  op <- prepare_conditional_local_weights(gaussian$V, i, N, label, gaussian$factor)
  op$intercept <- gaussian$mean[i] - sum(op$weights * gaussian$mean[N])
  op$coefficient <- as.numeric(gaussian$M[i, ] - colSums(op$weights * gaussian$M[N, , drop = FALSE]))
  op
}
