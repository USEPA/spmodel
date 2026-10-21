# iid models do not need spatial conditioning
get_conditional_iid_variance <- function(object, simulate_covparams = FALSE) {
  if (isTRUE(simulate_covparams) ||
      !is.null(object$random) || !inherits(object, c("splm", "spglm"))) return(NULL)
  spcov <- coef(object, type = "spcov")
  if (!inherits(spcov, c("none", "ie"))) return(NULL)
  variance <- spcov_ie_stabilized(spcov, object$diagtol)
  if (!is.finite(variance) || variance < 0) stop("Invalid ie variance.", call. = FALSE)
  if (variance == 0) return(NULL)
  variance
}

draw_conditional_iid_error <- function(n, samples, variance, chunk_size = 1000L) {
  value <- matrix(NA_real_, n, samples)
  for (start in seq.int(1L, samples, by = chunk_size)) {
    cols <- seq.int(start, min(samples, start + chunk_size - 1L))
    value[, cols] <- sqrt(variance) * matrix(rnorm(n * length(cols)), n, length(cols))
  }
  value
}

get_conditional_iid_glm <- function(object, Xnew, variance, samples, latent = FALSE, chunk_size = 1000L, center = get_latent_center(object)) {
  b <- center$betahat
  delta <- t(chol(vcov(object))) %*% matrix(rnorm(length(b) * samples), length(b), samples)
  beta <- sweep(delta, 1L, b, "+")
  rownames(beta) <- names(b)
  observed <- NULL
  if (latent) {
    X <- model.matrix(object)
    eta <- fitted(object, type = "link")
    D <- diag(get_D(object$family, eta, object$y, object$size,
      as.vector(coef(object, type = "dispersion"))))
    precision <- 1 / variance - D
    if (any(!is.finite(precision) | precision <= 0)) {
      stop("The conditional GLM latent covariance is not positive definite; check the fitted model and its convergence.", call. = FALSE)
    }
    q <- 1 / precision
    A <- (q / variance) * X
    w <- w_offset_free(center$eta, model.offset(model.frame(object)))
    observed <- matrix(NA_real_, NROW(X), samples)
    for (start in seq.int(1L, samples, by = chunk_size)) {
      cols <- seq.int(start, min(samples, start + chunk_size - 1L))
      observed[, cols] <- sweep(A %*% delta[, cols, drop = FALSE] +
        sqrt(q) * matrix(rnorm(NROW(X) * length(cols)), NROW(X), length(cols)), 1L, w, "+")
    }
  }
  value <- draw_conditional_iid_error(NROW(Xnew), samples, variance, chunk_size)
  for (start in seq.int(1L, samples, by = chunk_size)) {
    cols <- seq.int(start, min(samples, start + chunk_size - 1L))
    value[, cols] <- Xnew %*% beta[, cols, drop = FALSE] + value[, cols, drop = FALSE]
  }
  rownames(value) <- rownames(Xnew)
  list(beta = beta, newdata = value, latent = observed)
}
