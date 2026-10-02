#' Prepare a GLM neighborhood conditional given coefficients and earlier draws
#' @param covariance Spatial covariance, observed rows followed by target and previous rows.
#' @param X Observed neighborhood design.
#' @param Xnew Target and previous prediction design, in that order.
#' @param w Offset-free fitted latent values at observed neighbors.
#' @param D Observed response Hessian diagonal at fitted link values.
#' @param betahat Fitted coefficient vector.
#' @return Fixed intercept, coefficient/previous-value weights, and variance.
#' @noRd
get_conditional_glm_neighborhood <- function(covariance, X, Xnew, w, D, betahat) {
  n <- NROW(X)
  target <- n + seq_len(NROW(Xnew))
  S <- covariance[target, target, drop = FALSE]
  mu <- as.numeric(Xnew %*% betahat)
  F <- Xnew
  if (n) {
    observed <- seq_len(n)
    cov_lowchol <- t(chol(covariance[observed, observed, drop = FALSE]))
    Z <- forwardsolve(cov_lowchol, X)
    mH <- diag(n) - crossprod(cov_lowchol, D * cov_lowchol)
    conditional_chol <- tryCatch(chol(mH), error = function(e) {
      stop("The neighborhood latent precision is not positive definite; increase local$size or check the fitted model.", call. = FALSE)
    })
    E <- forwardsolve(cov_lowchol, covariance[observed, target, drop = FALSE])
    M <- backsolve(conditional_chol, forwardsolve(t(conditional_chol), Z))
    uncertainty <- forwardsolve(t(conditional_chol), E)
    S <- S - crossprod(E) + crossprod(uncertainty)
    mu <- mu + as.numeric(crossprod(E, forwardsolve(cov_lowchol, w) - Z %*% betahat))
    F <- Xnew - crossprod(E, Z) + crossprod(E, M)
  }
  previous <- seq_len(NROW(Xnew))[-1L]
  weights <- numeric()
  variance <- S[1L, 1L]
  intercept <- mu[1L]
  coefficient <- F[1L, ]
  if (length(previous)) {
    previous_chol <- chol(S[previous, previous, drop = FALSE])
    weights <- as.numeric(backsolve(previous_chol,
      forwardsolve(t(previous_chol), S[previous, 1L, drop = FALSE])))
    intercept <- intercept - sum(weights * mu[previous])
    coefficient <- coefficient - as.numeric(crossprod(weights, F[previous, , drop = FALSE]))
    variance <- variance - sum(weights * S[previous, 1L])
  }
  if (variance < -sqrt(.Machine$double.eps) * max(abs(S[1L, 1L]), .Machine$double.eps)) {
    stop("The neighborhood conditional variance is negative; check the fitted model.", call. = FALSE)
  }
  list(intercept = intercept, coefficient = as.numeric(coefficient),
    weights = weights, variance = max(variance, 0))
}

#' Prepare fixed Vecchia GLM operators without full observed latent draws
#' @param object A fitted spglm object.
#' @param newdata Processed prediction data in original row order.
#' @param Xnew Aligned prediction design matrix.
#' @param local_list Resolved Vecchia settings.
#' @return Ordered neighborhood operators and shared coefficient factor.
#' @noRd
get_conditional_vecchia_glm_context <- function(object, newdata, Xnew, local_list) {
  n <- NROW(object$obdata)
  m <- NROW(newdata)
  ord <- local_list$order
  X <- model.matrix(object)
  betahat <- coef(object)
  coefficient_factor <- t(chol(vcov(object)))
  eta <- fitted(object, type = "link")
  w <- w_offset_free(eta, model.offset(model.frame(object)))
  D <- diag(get_D(object$family, eta, object$y, object$size,
    as.vector(coef(object, type = "dispersion"))))
  columns <- unique(c(object$xcoord, object$ycoord, all.vars(object$random), all.vars(object$partition_factor)))
  pool <- rbind(object$obdata[, columns, drop = FALSE], newdata[ord, columns, drop = FALSE])
  covariance <- object[c("coefficients", "random", "partition_factor", "anisotropy",
    "xcoord", "ycoord", "dim_coords", "diagtol")]
  class(covariance) <- class(object)
  coords <- pool[, c(object$xcoord, object$ycoord), drop = FALSE]
  if (object$anisotropy) {
    spcov <- coef(object, type = "spcov")
    transformed <- transform_anis(pool, object$xcoord, object$ycoord, spcov[["rotate"]], spcov[["scale"]])
    coords <- data.frame(x = transformed$xcoord_val, y = transformed$ycoord_val)
  }
  operators <- vector("list", m)
  for (i in seq_len(m)) {
    candidates <- seq_len(n + i - 1L)
    if (local_list$method != "all" && length(candidates) > local_list$size) {
      if (local_list$method == "distance") {
        distance <- as.numeric(spdist_vectors2(coords[n + i, 1L], coords[n + i, 2L],
          coords[candidates, 1L], coords[candidates, 2L], sparse = FALSE))
        neighbors <- order(distance)[seq_len(local_list$size)]
      } else {
        covariance$obdata <- pool[candidates, , drop = FALSE]
        association <- covmatrix(covariance, pool[n + i, , drop = FALSE], cov_type = "pred.obs")
        neighbors <- order(abs(as.numeric(association)), decreasing = TRUE)[seq_len(local_list$size)]
      }
    } else {
      neighbors <- candidates
    }
    observed <- neighbors[neighbors <= n]
    previous <- neighbors[neighbors > n] - n
    prediction <- c(i, previous)
    covariance$obdata <- pool[c(observed, n + prediction), , drop = FALSE]
    operator <- get_conditional_glm_neighborhood(covmatrix(covariance),
      X[observed, , drop = FALSE], Xnew[ord[prediction], , drop = FALSE],
      w[observed], D[observed], betahat)
    operator$previous <- previous
    operators[[i]] <- operator
  }
  list(operators = operators, order = ord, betahat = betahat, coefficient_factor = coefficient_factor)
}

#' Simulate Vecchia GLM latent values using cached neighborhood conditionals
#' @param object A fitted spglm object.
#' @param newdata Processed prediction data in original row order.
#' @param Xnew Aligned prediction design matrix.
#' @param local_list Resolved Vecchia settings.
#' @param samples Number of draws.
#' @return Coefficient and offset-free prediction matrices.
#' @noRd
get_conditional_vecchia_glm <- function(object, newdata, Xnew, local_list, samples) {
  chunk_size <- get_conditional_glm_chunk_size(local_list)
  context <- get_conditional_vecchia_glm_context(object, newdata, Xnew, local_list)
  m <- NROW(newdata)
  p <- length(context$betahat)
  beta <- matrix(NA_real_, p, samples, dimnames = list(names(context$betahat), NULL))
  new_val <- matrix(NA_real_, m, samples, dimnames = list(rownames(Xnew), NULL))
  for (start in seq.int(1L, samples, by = chunk_size)) {
    columns <- seq.int(start, min(samples, start + chunk_size - 1L))
    count <- length(columns)
    delta <- context$coefficient_factor %*% matrix(rnorm(p * count), p, count)
    beta[, columns] <- sweep(delta, 1L, context$betahat, "+")
    draws <- matrix(NA_real_, m, count)
    for (i in seq_len(m)) {
      operator <- context$operators[[i]]
      mu <- operator$intercept + as.numeric(crossprod(operator$coefficient, delta))
      if (length(operator$previous)) {
        mu <- mu + as.numeric(crossprod(operator$weights, draws[operator$previous, , drop = FALSE]))
      }
      draws[i, ] <- mu + sqrt(operator$variance) * rnorm(count)
    }
    new_val[context$order, columns] <- draws
  }
  list(beta = beta, newdata = new_val)
}
