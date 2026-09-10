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
  SigInv <- Matrix::bdiag(SigInv_list)
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
  rt_SigInv_r <- crossprod(r, SigInv) %*% r

  # get wolfinger objects
  y <- as.vector(do.call("rbind", data_object$y_list))
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }
  # l00 is minus twice the conditional data log-likelihood at the converged w,
  # l01 is the Laplace correction (log determinant of the negative Hessian at
  # the mode); l1/l2/(l3) extend the usual Gaussian wolfinger pieces so that
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
  # l01 is the Laplace correction (log determinant of the negative Hessian at the mode)
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
#'   shape), \code{mHldet} (the log-determinant of the negative Hessian), and,
#'   if \code{ret_mHInv} is \code{TRUE}, \code{mHInv} (the inverse of the
#'   negative Hessian). When there is more than one partition, the update
#'   uses the Sherman-Morrison-Woodbury identity (via \code{smw_HInv()})
#'   instead of a direct solve, since the Hessian is otherwise too large to invert
#'
#' @noRd
get_w_and_H_spglm <- function(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
  family <- data_object$family
  SigInv <- Matrix::bdiag(SigInv_list)
  # Ptheta is the precision matrix projected off the fixed-effect space
  # (SigInv adjusted for estimating betahat), used in the score/Hessian below
  Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
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

    mHldet <- as.numeric(determinant(-H, logarithm = TRUE)$modulus)
    w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
    if (ret_mHInv) {
      # not done above because this is only for model stats and solve(H) slower than solve(H, g)
      HInv <- solve(H)
      w_and_H_list$mHInv <- -HInv
    }
  } else {
    # multi-partition case: the full Hessian is too large to invert directly,
    # so its block-diagonal-plus-low-rank structure is exploited via the
    # Sherman-Morrison-Woodbury identity instead (see smw_HInv()/smw_mHldet())
    # add cov_betahat_Inv stability by same diagonal tolerance as this can have problems too
    diag(cov_betahat_Inv) <- diag(cov_betahat_Inv) + data_object$diagtol

    while (iter < 50 && wdiffmax > 1e-4) {
      iter <- iter + 1
      # compute the d vector
      d <- get_d(family, w + offset, y, size, dispersion)
      # and then the gradient vector
      g <- d - Ptheta %*% w
      # Next, compute H
      D <- get_D(family, w + offset, y, size, dispersion)
      D_diag <- diag(D)
      # split the diagonal Hessian contribution back out by partition so each
      # partition's block of -H (D - SigInv) can be combined with that
      # partition's SigInv block below
      D_list <- lapply(split(D_diag, sort(data_object$local_index)), function(x) Diagonal(x = x))
      # cholesky products
      if (data_object$parallel) {
        cluster_list <- lapply(seq_along(D_list), function(l) {
          cluster_list_element <- list(
            D = D_list[[l]],
            S = SigInv_list[[l]]
          )
        })
        DSigInv_list <- parallel::parLapply(data_object$cl, cluster_list, get_DSigInv_parallel)
        names(DSigInv_list) <- names(D_list)
      } else {
        DSigInv_list <- mapply(
          D = D_list, S = SigInv_list,
          function(D, S) get_DSigInv(D, S),
          SIMPLIFY = FALSE
        )
      }

      if (data_object$parallel) {
        cluster_list <- DSigInv_list
        DSigInv_Inv_list <- parallel::parLapply(data_object$cl, cluster_list, solve)
        names(DSigInv_Inv_list) <- names(D_list)
      } else {
        DSigInv_Inv_list <- lapply(DSigInv_list, function(x) solve(x))
      }
      DSigInv_Inv <- Matrix::bdiag(DSigInv_Inv_list)
      HInv <- smw_HInv(AInv = DSigInv_Inv, U = SigInv_X, CInv = cov_betahat_Inv)
      solveHg <- HInv %*% g
      wnew <- w - solveHg
      # check overshoot on loglik surface
      dnew <- get_d(family, wnew + offset, y, size, dispersion)
      gnew <- dnew - Ptheta %*% wnew
      if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing ie at a known, non-zero value (via spcov_initial), or fixing dispersion at one (via dispersion_initial).", call. = FALSE)
      if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
      wdiffmax <- max(abs(wnew - w))
      # update w
      w <- wnew
    }

    mHldet <- smw_mHldet(A_list = DSigInv_list, AInv = DSigInv_Inv, U = SigInv_X, C = cov_betahat, CInv = cov_betahat_Inv)
    w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
    if (ret_mHInv) {
      w_and_H_list$mHInv <- -HInv
    }
  }

  w_and_H_list
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
#'   shape), \code{mHldet} (the log-determinant of the negative Hessian), and,
#'   if \code{ret_mHInv} is \code{TRUE}, \code{mHInv} (the inverse of the
#'   negative Hessian)
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

  mHldet <- as.numeric(determinant(-H, logarithm = TRUE)$modulus)
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

#' Invert the negative Hessian via the Sherman-Morrison-Woodbury identity
#'
#' @param AInv The inverse of the block-diagonal part of \eqn{-H}
#' @param U The (tall) coupling matrix (\code{SigInv_X})
#' @param CInv The inverse of the low-rank part (\code{cov_betahat_Inv})
#'
#' @return The inverse of \eqn{-H = AInv^{-1} - U C U'}, computed without
#'   forming or inverting the full (partition-sized) matrix directly
#'
#' @noRd
smw_HInv <- function(AInv, U, CInv) {
  # "mid" is the small (p x p, p = number of fixed effects) matrix that has
  # to be inverted, in place of inverting the full (n x n) -H
  mid <- CInv + t(U) %*% AInv %*% U
  AInv - (AInv %*% U) %*% solve(mid) %*% (t(U) %*% AInv)
}

#' Log-determinant of the negative Hessian via the matrix determinant lemma
#'
#' @param A_list A list of the block-diagonal parts of \eqn{-H}, one per partition
#' @param AInv The inverse of the block-diagonal part of \eqn{-H}
#' @param U The (tall) coupling matrix (\code{SigInv_X})
#' @param C The low-rank part (\code{cov_betahat})
#' @param CInv The inverse of the low-rank part (\code{cov_betahat_Inv})
#'
#' @return The log-determinant of \eqn{-H}, computed without forming or
#'   taking the determinant of the full (partition-sized) matrix directly
#'
#' @noRd
smw_mHldet <- function(A_list, AInv, U, C, CInv) {
  Aldet <- sum(unlist(lapply(A_list, function(x) determinant(x, logarithm = TRUE)$modulus))) # must be positive det for -H
  Cldet <- 2 * sum(log(diag(t(chol(C)))))
  mid <- CInv + t(U) %*% AInv %*% U
  midldet <- determinant(mid, logarithm = TRUE)$modulus
  as.numeric(Aldet + Cldet + midldet)
}

#' Compute the block-diagonal part of the negative Hessian for one partition
#'
#' @param D The (diagonal) GLM Hessian contribution for the partition
#' @param SigInv The partition's inverse covariance matrix
#'
#' @return \code{D - SigInv}, the partition's contribution to \eqn{-H}
#'   (before the low-rank \code{U C U'} correction)
#'
#' @noRd
get_DSigInv <- function(D, SigInv) {
  D - SigInv
}

#' Parallel-friendly wrapper around \code{get_DSigInv()}
#'
#' @param cluster_list A list with elements \code{D} and \code{S}
#'
#' @return The same value as \code{get_DSigInv()}, for use with \code{parallel::parLapply()}
#'
#' @noRd
get_DSigInv_parallel <- function(cluster_list) {
  D <- cluster_list$D
  S <- cluster_list$S
  get_DSigInv(D, S)
}
