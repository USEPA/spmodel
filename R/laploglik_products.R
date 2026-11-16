#' Find relevant products to use in Gaussian log-likelihood calculations
#'
#' @param spcov_params_val A \code{spcov_params} object
#' @param dispersion_params_val A \code{dispersion_params} object
#' @param ... other arguments
#'
#' @return The relevant Gaussian log-likelihood products
#'
#' @noRd
# dispatches on the covariance function class, mirroring gloglik_products()
# but for GLM responses: the latent Gaussian random effect w is integrated
# out with a Laplace approximation rather than observed directly
laploglik_products <- function(spcov_params_val, dispersion_params_val, ...) {
  UseMethod("laploglik_products", spcov_params_val)
}
#' @export
laploglik_products.exponential <- function(spcov_params_val, dispersion_params_val, data_object, estmethod,
                                           dist_matrix_list, randcov_params_val, ...) {
  # making a covariance matrix
  cov_matrix_list <- get_cov_matrix_list(spcov_params_val, dist_matrix_list, randcov_params_val, data_object$randcov_list, data_object$partition_list,
    diagtol = data_object$diagtol
  )


  # cholesky products
  # cov_matrix_list holds one block per big-data partition (a single block
  # when there is no partitioning); each block's Cholesky factorization is
  # independent, so it is parallelized across a cluster when requested
  if (data_object$parallel) {
    cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
      cluster_list_element <- list(
        c = cov_matrix_list[[l]],
        x = data_object$X_list[[l]],
        y = data_object$y_list[[l]]
      )
    })
    cholprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_cholprods_glm_parallel)
    names(cholprods_list) <- names(cov_matrix_list)
  } else {
    cholprods_list <- mapply(
      c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
      function(c, x, y) get_cholprods_glm(c, x, y),
      SIMPLIFY = FALSE
    )
  }

  SigInv_list <- lapply(cholprods_list, function(x) x$SigInv)
  SigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SigInv_X))


  # storing relevant products
  ## lower chol %*% X
  SqrtSigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_X))
  ## lower chol %*% y
  SqrtSigInv_y <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_y))
  # covariance of beta hat
  ## t(X) %*% sigma_inverse %*% X
  Xt_SigInv_X <- crossprod(SqrtSigInv_X, SqrtSigInv_X)
  ## t(X) %*% sigma_inverse %*% X)^(-1)
  Xt_SigInv_X_upchol <- chol(Xt_SigInv_X)
  cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

  # find dispersion
  dispersion <- as.vector(dispersion_params_val) # take class away

  # newton rhapson
  # find the latent w that maximizes the joint (data + random effect)
  # log-likelihood -- this is the mode used by the Laplace approximation to
  # the marginal (w integrated out) likelihood
  w_and_H <- get_w_and_H_spglm(
    data_object, dispersion,
    SigInv_list, SigInv_X, cov_betahat, Xt_SigInv_X, estmethod
  )

  w <- w_and_H$w
  mHldet <- w_and_H$mHldet

  betahat <- tcrossprod(cov_betahat, SigInv_X) %*% w
  X <- do.call("rbind", data_object$X_list)
  r <- w - X %*% betahat
  sizes <- vapply(SigInv_list, nrow, 1L)
  idx <- split(seq_len(sum(sizes)), rep(seq_along(sizes), sizes))
  rt_SigInv_r <- sum(mapply(function(S, i) as.numeric(crossprod(r[i], S %*% r[i])), SigInv_list, idx))

  # get wolfinger objects
  y <- as.vector(do.call("rbind", data_object$y_list))
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }
  # l00 is minus twice the conditional data log-likelihood at the converged w,
  # l01 is the Laplace log determinant; l1/l2/(l3) extend the Gaussian pieces so that
  # get_minustwolaploglik() can combine them the same way as get_minustwologlik()
  l00 <- get_l00(data_object$family, w, y, data_object$size, dispersion)
  l01 <- mHldet
  l1 <- sum(unlist(lapply(cholprods_list, function(x) 2 * sum(log(diag(x$Sig_lowchol))))))
  l2 <- as.numeric(rt_SigInv_r)

  # returning relevant quantities
  if (estmethod == "reml") {
    l3 <- 2 * sum(log(diag(Xt_SigInv_X_upchol)))
    return(list(l00 = l00, l01 = l01, l1 = l1, l2 = l2, l3 = l3))
  }

  if (estmethod == "ml") {
    return(list(l00 = l00, l01 = l01, l1 = l1, l2 = l2))
  }
}
#' @export
laploglik_products.spherical <- laploglik_products.exponential
#' @export
laploglik_products.gaussian <- laploglik_products.exponential
#' @export
laploglik_products.triangular <- laploglik_products.exponential
#' @export
laploglik_products.circular <- laploglik_products.exponential
#' @export
laploglik_products.none <- laploglik_products.exponential
#' @export
laploglik_products.ie <- laploglik_products.none
#' @export
laploglik_products.cubic <- laploglik_products.exponential
#' @export
laploglik_products.pentaspherical <- laploglik_products.exponential
#' @export
laploglik_products.cosine <- laploglik_products.exponential
#' @export
laploglik_products.wave <- laploglik_products.exponential
#' @export
laploglik_products.jbessel <- laploglik_products.exponential
#' @export
laploglik_products.gravity <- laploglik_products.exponential
#' @export
laploglik_products.rquad <- laploglik_products.exponential
#' @export
laploglik_products.magnetic <- laploglik_products.exponential

#' @export
laploglik_products.matern <- laploglik_products.exponential
#' @export
laploglik_products.cauchy <- laploglik_products.exponential
#' @export
laploglik_products.pexponential <- laploglik_products.exponential

#' @export
laploglik_products.car <- function(spcov_params_val, dispersion_params_val, data_object, estmethod,
                                   dist_matrix_list, randcov_params_val, ...) {
  # car/sar models parameterize the *precision* (inverse covariance) matrix
  # directly and sparsely, so SigInv and its log determinant come from a
  # dedicated helper rather than from Cholesky-factoring a dense Sigma
  spautor_cov_matrixInv_val <- spautor_cov_matrixInv(
    spcov_params_val, data_object,
    dist_matrix_list, randcov_params_val
  )

  SigInv <- spautor_cov_matrixInv_val$SigInv
  Sigldet <- spautor_cov_matrixInv_val$Sigldet

  # finding relevant quantities for likelihood
  SigInv_X <- SigInv %*% data_object$X
  Xt_SigInv_X <- crossprod(data_object$X, SigInv_X)
  Xt_SigInv_X_upchol <- chol(forceSymmetric(Xt_SigInv_X))
  cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

  # find dispersion
  dispersion <- as.vector(dispersion_params_val) # take class away

  # newton rhapson
  # find the latent w that maximizes the joint (data + random effect)
  # log-likelihood -- this is the mode used by the Laplace approximation to
  # the marginal (w integrated out) likelihood
  w_and_H <- get_w_and_H_spgautor(
    data_object, dispersion,
    SigInv, SigInv_X, cov_betahat, Xt_SigInv_X, estmethod
  )

  w <- w_and_H$w
  mHldet <- w_and_H$mHldet

  betahat <- tcrossprod(cov_betahat, SigInv_X) %*% w

  X <- data_object$X
  r <- w - X %*% betahat
  rt_SigInv_r <- crossprod(r, SigInv) %*% r

  # get wolfinger objects
  y <- data_object$y
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }
  # l00 is minus twice the conditional data log-likelihood at the converged w,
  # l01 is the Laplace log determinant
  l00 <- get_l00(data_object$family, w, y, data_object$size, dispersion)
  l01 <- mHldet
  l1 <- Sigldet
  l2 <- as.numeric(rt_SigInv_r)


  # returning relevant quantities
  if (estmethod == "reml") {
    l3 <- 2 * sum(log(diag(Xt_SigInv_X_upchol)))
    return(list(l00 = l00, l01 = l01, l1 = l1, l2 = l2, l3 = l3))
  }

  if (estmethod == "ml") {
    return(list(l00 = l00, l01 = l01, l1 = l1, l2 = l2))
  }
}
#' @export
laploglik_products.sar <- laploglik_products.car


#' Newton-Raphson solve for the latent \code{w} vector (\code{spglm()} models)
#'
#' @param data_object The data object
#' @param dispersion The dispersion parameter
#' @param SigInv_list A list of partition-wise inverse covariance matrices
#' @param SigInv_X \code{SigInv \%*\% X}
#' @param cov_betahat The covariance matrix of betahat
#' @param cov_betahat_Inv The inverse of \code{cov_betahat} (i.e. \eqn{X'\Sigma^{-1}X})
#' @param estmethod The estimation method
#' @param ret_mHInv Whether to also return the inverse of the negative Hessian
#'
#' @return A list with elements \code{w} (the converged latent predictor
#'   vector), \code{H} (always \code{NULL}; retained for a consistent return
#'   shape), \code{mHldet} (the Laplace log determinant, unprojected under ML
#'   and projected under REML), and,
#'   if \code{ret_mHInv} is \code{TRUE}, \code{mHInv} (the inverse of the
#'   projected negative Hessian, represented by block factors for multiple partitions).
#'   When there is more than one partition, the update
#'   applies the Sherman-Morrison-Woodbury identity through block solves
#'   instead of forming the full inverse in every Newton iteration
#'
#' @noRd
get_w_and_H_spglm <- function(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
  family <- data_object$family
  y <- as.vector(do.call("rbind", data_object$y_list))
  size <- data_object$size
  w <- get_w_init(family, y, dispersion)
  wdiffmax <- Inf
  iter <- 0

  # The offset is a known, non-estimated shift on the link scale, so w and the
  # linear predictor are not the same vector: w is the offset-free latent
  # process the optimizer solves for (and the scale of the
  # spatial covariance), while the family log-likelihood is always evaluated at
  # the linear predictor w + offset. Every get_d()/get_D() call below therefore
  # takes w + offset, and every Ptheta product takes w alone. Defaulting the
  # offset to 0 keeps that distinction visible in one place instead of
  # duplicating it across if/else branches, matching glm().
  offset <- if (is.null(data_object$offset)) 0 else as.vector(data_object$offset)

  # single-partition case: the Hessian is small enough to solve directly
  if (length(SigInv_list) == 1) {
    SigInv <- as.matrix(SigInv_list[[1]])
    SigInv_X <- as.matrix(SigInv_X)
    cov_betahat <- as.matrix(cov_betahat)
    Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
    while (iter < 50 && wdiffmax > 1e-4) {
      iter <- iter + 1
      # compute the d vector
      d <- get_d(family, w + offset, y, size, dispersion)
      # and then the gradient vector
      g <- as.vector(d - Ptheta %*% w)
      # Next, compute H
      D <- get_D(family, w + offset, y, size, dispersion)
      mH <- Ptheta
      diag(mH) <- diag(mH) - as.vector(diag(D))
      R <- tryCatch(chol(mH), error = function(e) NULL)
      solveHg <- if (is.null(R)) solve(-mH, g) else -backsolve(R, forwardsolve(t(R), g))
      solveHg <- as.vector(solveHg)
      wnew <- w - solveHg
      # check overshoot on loglik surface
      dnew <- get_d(family, wnew + offset, y, size, dispersion)
      gnew <- dnew - Ptheta %*% wnew
      if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing ie at a known, non-zero value (via spcov_initial), or fixing dispersion at one (via dispersion_initial).", call. = FALSE)
      if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
      wdiffmax <- max(abs(wnew - w))
      w <- wnew
    }

    # ML holds beta fixed in the integral; Newton uses the profiled Hessian.
    if (estmethod == "ml") {
      m <- SigInv
      diag(m) <- diag(m) - as.vector(diag(D))
      Rm <- tryCatch(chol(m), error = function(e) NULL)
      mHldet <- if (is.null(Rm)) as.numeric(determinant(m, logarithm = TRUE)$modulus) else 2 * sum(log(diag(Rm)))
    } else {
      mHldet <- if (is.null(R)) as.numeric(determinant(mH, logarithm = TRUE)$modulus) else 2 * sum(log(diag(R)))
    }
    w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
    if (ret_mHInv) {
      w_and_H_list$mHInv <- if (is.null(R)) solve(mH) else chol2inv(R)
    }
  } else {
    SigInv_list <- lapply(SigInv_list, as.matrix)
    SigInv_X <- as.matrix(SigInv_X)
    cov_betahat <- as.matrix(cov_betahat)
    cov_betahat_Inv <- as.matrix(cov_betahat_Inv)
    sizes <- vapply(SigInv_list, nrow, 1L)
    idx <- split(seq_len(sum(sizes)), rep(seq_along(sizes), sizes))
    stopifnot(identical(unname(lengths(idx)),
      unname(lengths(split(seq_along(y), sort(data_object$local_index))))))
    Pw <- function(v) {
      Sv <- unlist(Map(function(S, i) S %*% v[i], SigInv_list, idx), use.names = FALSE)
      Sv - as.vector(SigInv_X %*% (cov_betahat %*% crossprod(SigInv_X, v)))
    }
    p <- ncol(SigInv_X)
    # add cov_betahat_Inv stability by same diagonal tolerance as this can have problems too
    diag(cov_betahat_Inv) <- diag(cov_betahat_Inv) + data_object$diagtol

    while (iter < 50 && wdiffmax > 1e-4) {
      iter <- iter + 1
      # compute the d vector
      d <- get_d(family, w + offset, y, size, dispersion)
      # and then the gradient vector
      g <- d - Pw(w)
      # Next, compute H
      D <- get_D(family, w + offset, y, size, dispersion)
      F <- factor_laplace_blocks(SigInv_list, as.vector(diag(D)), idx,
        data_object$parallel, data_object$cl)
      AUg <- apply_Ainv(F, cbind(SigInv_X, g), idx)
      AU <- AUg[, seq_len(p), drop = FALSE]
      Ag <- AUg[, p + 1L]
      mid <- cov_betahat_Inv + crossprod(SigInv_X, AU)
      # Woodbury product; see technical.Rmd #sec:bigdata-glm.
      solveHg <- as.vector(Ag - AU %*% solve(mid, crossprod(SigInv_X, Ag)))
      wnew <- w - solveHg
      # check overshoot on loglik surface
      dnew <- get_d(family, wnew + offset, y, size, dispersion)
      gnew <- dnew - Pw(wnew)
      if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing ie at a known, non-zero value (via spcov_initial), or fixing dispersion at one (via dispersion_initial).", call. = FALSE)
      if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
      wdiffmax <- max(abs(wnew - w))
      # update w
      w <- wnew
    }

    mHldet <- logdet_A(F)
    if (estmethod == "reml") {
      mHldet <- mHldet + 2 * sum(log(diag(chol(cov_betahat)))) +
        as.numeric(determinant(mid, logarithm = TRUE)$modulus)
    }
    w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
    if (ret_mHInv) {
      w_and_H_list$mHInv <- structure(list(F = F, idx = idx, AU = AU, mid_inv = solve(mid)),
        class = "laplace_mhinv")
    }
  }

  w_and_H_list
}

factor_laplace_block <- function(block) {
  m <- block$S
  diag(m) <- diag(m) - block$D
  R <- tryCatch(chol(m), error = function(e) NULL)
  if (is.null(R)) list(type = "lu", A = -m) else list(type = "chol", R = R)
}

factor_laplace_blocks <- function(S_list, D, idx, parallel, cl) {
  blocks <- Map(function(S, i) list(S = S, D = D[i]), S_list, idx)
  if (parallel) parallel::parLapply(cl, blocks, factor_laplace_block) else lapply(blocks, factor_laplace_block)
}

apply_Ainv <- function(F, V, idx) {
  V <- as.matrix(V)
  do.call(rbind, Map(function(f, i) {
    if (f$type == "chol") {
      -backsolve(f$R, forwardsolve(t(f$R), V[i, , drop = FALSE]))
    } else {
      solve(f$A, V[i, , drop = FALSE])
    }
  }, F, idx))
}

logdet_A <- function(F) {
  sum(vapply(F, function(f) {
    if (f$type == "chol") 2 * sum(log(diag(f$R))) else as.numeric(determinant(f$A, logarithm = TRUE)$modulus)
  }, numeric(1)))
}

mhinv_apply <- function(M, V) {
  if (!inherits(M, "laplace_mhinv")) return(as.matrix(M %*% V))
  V <- as.matrix(V)
  -(apply_Ainv(M$F, V, M$idx) - M$AU %*% (M$mid_inv %*% crossprod(M$AU, V)))
}

mhinv_diag <- function(M) {
  if (!inherits(M, "laplace_mhinv")) return(Matrix::diag(M))
  diag_Ainv <- unlist(lapply(M$F, function(f) {
    if (f$type == "chol") -diag(chol2inv(f$R)) else diag(solve(f$A))
  }), use.names = FALSE)
  -(diag_Ainv - rowSums((M$AU %*% M$mid_inv) * M$AU))
}

#' Newton-Raphson solve for the latent \code{w} vector (\code{spgautor()} models)
#'
#' @param data_object The data object
#' @param dispersion The dispersion parameter
#' @param SigInv The inverse covariance matrix
#' @param SigInv_X \code{SigInv \%*\% X}
#' @param cov_betahat The covariance matrix of betahat
#' @param cov_betahat_Inv The inverse of \code{cov_betahat} (i.e. \eqn{X'\Sigma^{-1}X})
#' @param estmethod The estimation method
#' @param ret_mHInv Whether to also return the inverse of the negative Hessian
#'
#' @return A list with elements \code{w} (the converged latent predictor
#'   vector), \code{H} (always \code{NULL}; retained for a consistent return
#'   shape), \code{mHldet} (the Laplace log determinant, unprojected under ML
#'   and projected under REML), and,
#'   if \code{ret_mHInv} is \code{TRUE}, \code{mHInv} (the inverse of the
#'   projected negative Hessian)
#'
#' @noRd
get_w_and_H_spgautor <- function(data_object, dispersion, SigInv, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
  family <- data_object$family
  Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
  y <- data_object$y
  size <- data_object$size
  w <- get_w_init(family, y, dispersion)
  wdiffmax <- Inf
  iter <- 0

  # see get_w_and_H_spglm(): w is the offset-free latent process, while the
  # family log-likelihood is evaluated at the linear predictor w + offset
  offset <- if (is.null(data_object$offset)) 0 else as.vector(data_object$offset)

  while (iter < 50 && wdiffmax > 1e-4) {
    iter <- iter + 1
    # compute the d vector
    d <- get_d(family, w + offset, y, size, dispersion)
    # and then the gradient vector
    g <- d - Ptheta %*% w
    # Next, compute H
    D <- get_D(family, w + offset, y, size, dispersion)
    H <- D - Ptheta # not PD but -H is
    solveHg <- solve(H, g)
    wnew <- w - solveHg
    # check overshoot on loglik surface
    dnew <- get_d(family, wnew + offset, y, size, dispersion)
    gnew <- dnew - Ptheta %*% wnew
    if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing ie at a known, non-zero value (via spcov_initial), or fixing dispersion at one (via dispersion_initial).", call. = FALSE)
    if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
    wdiffmax <- max(abs(wnew - w))
    w <- wnew
  }

  mHldet <- as.numeric(determinant(if (estmethod == "ml") SigInv - D else -H, logarithm = TRUE)$modulus)
  w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
  if (ret_mHInv) {
    # not done above because this is only for model stats and solve(H) slower than solve(H, g)
    HInv <- solve(H)
    w_and_H_list$mHInv <- -HInv
  }

  w_and_H_list
}

#' Compute the gradient of the Laplace log-likelihood with respect to \code{w}
#'
#' @param family The response family
#' @param w The latent (link-scale) predictor vector
#' @param y Response vector
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter
#'
#' @return The gradient vector (denoted \eqn{d} in the package's Laplace
#'   approximation derivation)
#'
#' @noRd
get_d <- function(family, w, y, size, dispersion) {
  if (family == "poisson") {
    d <- -exp(w) + y
  } else if (family == "nbinomial") {
    d <- dispersion * (y - exp(w)) / (dispersion + exp(w))
  } else if (family == "binomial") {
    d <- y - size * expit(w)
  } else if (family == "Gamma") {
    d <- -dispersion + dispersion * y * exp(-w)
  } else if (family == "inverse.gaussian") {
    # d <- 1 / dispersion * (y - exp(w)) / exp(2 * w)
    d <- dispersion * (y / (2 * exp(w)) - exp(w) / (2 * y)) + 1 / 2
  } else if (family == "beta") {
    one_expw <- 1 + exp(w)
    k0 <- digamma(dispersion * exp(w) / one_expw) - digamma(dispersion / one_expw) + log(1 / y - 1)
    d <- -dispersion * exp(w) * k0 / one_expw^2
  }
  d
}

#' Compute the (diagonal) Hessian of the Laplace log-likelihood with respect to \code{w}
#'
#' @param family The response family
#' @param w The latent (link-scale) predictor vector
#' @param y Response vector
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter
#'
#' @return A diagonal matrix (denoted \eqn{D} in the package's Laplace
#'   approximation derivation), diagonal because observations are
#'   conditionally independent given \code{w}
#'
#' @noRd
get_D <- function(family, w, y, size, dispersion) {
  w <- as.vector(w)

  if (family == "poisson") {
    D_vec <- -exp(w)
  } else if (family == "nbinomial") {
    D_vec <- -(dispersion * exp(w) * (dispersion + y)) / ((dispersion + exp(w))^2)
  } else if (family == "binomial") {
    D_vec <- -size * expit(w) / (1 + exp(w))
  } else if (family == "Gamma") {
    D_vec <- -dispersion * y * exp(-w)
  } else if (family == "inverse.gaussian") {
    # D_vec <- 1 / dispersion * (exp(w) - 2 * y) / exp(2 * w)
    D_vec <- -dispersion * (exp(2 * w) + y^2) / (2 * y * exp(w))
  } else if (family == "beta") {
    one_expw <- 1 + exp(w)
    k0 <- digamma(dispersion * exp(w) / one_expw) - digamma(dispersion / one_expw) + log(1 / y - 1)
    # get_d() has d = -A(w) * k0 with A(w) = dispersion * exp(w) / one_expw^2,
    # so D = -A'(w) * k0 - A(w) * dk0/dw, which repackages into the -2 sinh(w)
    # k0 term plus the trigamma term below. k0 already ends in log((1 - y) / y). A previous version of this line added
    # 2 * atanh(1 - 2 * y) (algebraically the same quantity) a second time,
    # double-counting this piece in the second derivative.
    k1 <- dispersion * (trigamma(dispersion * exp(w) / one_expw) + trigamma(dispersion / one_expw)) - 2 * sinh(w) * k0
    D_vec <- -dispersion * exp(2 * w) * k1 / one_expw^4
  }
  D <- Diagonal(x = D_vec)
}

get_D3 <- function(family, w, y, size, dispersion) {
  w <- as.numeric(w)
  mu <- exp(w)
  if (family == "poisson") {
    value <- -mu
  } else if (family == "nbinomial") {
    p <- stats::plogis(w - log(dispersion))
    value <- -(dispersion + y) * p * (1 - p) * (1 - 2 * p)
  } else if (family == "binomial") {
    p <- expit(w)
    value <- -size * p * (1 - p) * (1 - 2 * p)
  } else if (family == "Gamma") {
    value <- dispersion * y * exp(-w)
  } else if (family == "inverse.gaussian") {
    value <- -dispersion * (mu - y^2 * exp(-w)) / (2 * y)
  } else if (family == "beta") {
    p <- expit(w)
    a <- p * (1 - p)
    b <- 1 - 2 * p
    k0 <- digamma(dispersion * p) - digamma(dispersion * (1 - p)) + log((1 - y) / y)
    k1 <- trigamma(dispersion * p) + trigamma(dispersion * (1 - p))
    k2 <- psigamma(dispersion * p, 2) - psigamma(dispersion * (1 - p), 2)
    value <- -dispersion * a * (b^2 - 2 * a) * k0 -
      3 * dispersion^2 * a^2 * b * k1 - dispersion^3 * a^3 * k2
  }
  as.numeric(value)
}

get_laplace_correction <- function(mHInv, family, eta, y, size, dispersion) {
  v <- mhinv_diag(mHInv)
  t3 <- get_D3(family, eta, y, size, dispersion)
  delta <- 0.5 * as.numeric(mhinv_apply(mHInv, t3 * v))
  list(delta = delta)
}

get_latent_center <- function(object, fit_correct = TRUE) {
  if (!is.logical(fit_correct) || length(fit_correct) != 1L || is.na(fit_correct)) {
    stop("fit_correct must be TRUE or FALSE.", call. = FALSE)
  }
  eta <- fitted(object, type = "link")
  active <- fit_correct && !is.null(object$laplace)
  if (!active) return(list(eta = eta, betahat = coef(object), active = FALSE))
  list(eta = eta + object$laplace$delta[names(eta)], betahat = object$laplace$betahat, active = TRUE)
}

#' Get a starting value for the Newton-Raphson solve of \code{w}
#'
#' @param family The response family
#' @param y Response vector
#' @param dispersion The dispersion parameter (unused; kept for a consistent signature)
#'
#' @return An initial guess for the latent (link-scale) predictor vector
#'
#' @noRd
get_w_init <- function(family, y, dispersion) {
  if (family == "poisson") {
    w_init <- 0.5 * log(y + 1)
  } else if (family == "nbinomial") {
    w_init <- 0.5 * log(y + 1)
  } else if (family == "binomial") {
    w_init <- rep(0, times = length(y))
  } else if (family == "Gamma") {
    w_init <- 0.5 * log(y + 1)
  } else if (family == "inverse.gaussian") {
    w_init <- 0.5 * log(y + 1)
  } else if (family == "beta") {
    w_init <- rep(0, times = length(y))
  }
  w_init
}

#' Compute minus twice the conditional log-likelihood \eqn{\log[y|g^{-1}(w),\phi]}
#'
#' @param family The response family
#' @param w The latent (link-scale) predictor vector
#' @param y Response vector
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter
#'
#' @return Minus twice the conditional (data-model) log-likelihood, evaluated
#'   at the converged \code{w}; one of the terms in the Laplace-approximated
#'   log-likelihood (denoted \eqn{l_{00}})
#'
#' @noRd
get_l00 <- function(family, w, y, size, dispersion) {
  w <- as.vector(w)
  y <- as.vector(y)
  # -2 is for -2ll constant
  if (family == "poisson") {
    mu <- exp(w)
    l00 <- -2 * sum(dpois(y, lambda = mu, log = TRUE))
  } else if (family == "nbinomial") {
    mu <- exp(w)
    l00 <- -2 * sum(dnbinom(x = y, mu = mu, size = dispersion, log = TRUE))
  } else if (family == "binomial") {
    mu <- expit(w)
    l00 <- -2 * sum(dbinom(y, size, mu, log = TRUE))
  } else if (family == "Gamma") {
    mu <- exp(w)
    # disp_recip <- 1 / dispersion
    # l00 <- -2 * sum(dgamma(y, shape = disp_recip, scale = dispersion * mu, log = TRUE))
    l00 <- -2 * sum(dgamma(y, shape = dispersion, scale = mu / dispersion, log = TRUE))
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    # matches statmod::dinvgauss(y, mean = mu, dispersion = 1 / (mu * dispersion), log = TRUE)
    # without requiring a statmod dependency
    l00 <- -2 * sum(1 / 2 * (log(dispersion) + log(exp(w)) - log(2 * pi) - log(y^3)) - dispersion * (y - exp(w))^2 / (2 * exp(w) * y))
  } else if (family == "beta") {
    mu <- expit(w)
    a <- mu * dispersion
    b <- (1 - mu) * dispersion
    l00 <- -2 * sum(dbeta(x = y, shape1 = a, shape2 = b, log = TRUE))
  }
  l00
}
