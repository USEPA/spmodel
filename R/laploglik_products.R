#' Find relevant products to use in Gaussian log-likelihood calculations
#'
#' @param spcov_params_val A \code{spcov_params} object
#' @param dispersion_params_val A \code{dispersion_params} object
#' @param ... other arguments
#'
#' @return The relevant Gaussian log-likelihood products
#'
#' @noRd
laploglik_products <- function(spcov_params_val, dispersion_params_val, ...) {
  UseMethod("laploglik_products", spcov_params_val)
}
#' @export
laploglik_products.exponential <- function(spcov_params_val, dispersion_params_val, data_object, estmethod,
                                         dist_matrix_list, randcov_params_val) {

  if (inherits(spcov_params_val, "none") && spcov_params_val[["ie"]] < 1e-4) {
    # instability when in smw of H when ie is small enough and the covariance is "none"
    # spcov_matrix.none does not adequately handle this when de = 0 because max set to 0
    # maybe another statement after the first max setting to check for zero to get
    # around this case?
    spcov_params_val[["ie"]] <- 1e-4
  }


  # making a covariance matrix
  cov_matrix_list <- get_cov_matrix_list(spcov_params_val, dist_matrix_list, randcov_params_val, data_object$randcov_list, data_object$partition_list)


  # cholesky products
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
  w_and_H <- get_w_and_H_spglm(data_object, dispersion,
                         SigInv_list, SigInv_X, cov_betahat, Xt_SigInv_X, estmethod)

  w <- w_and_H$w
  # H <- w_and_H$H
  mHldet <- w_and_H$mHldet

  betahat <- tcrossprod(cov_betahat, SigInv_X) %*% w
  X <- do.call("rbind", data_object$X_list)
  r <- w - X %*% betahat
  rt_SigInv_r <- crossprod(r, SigInv) %*% r

  # get wolfinger objects
  y <- as.vector(do.call("rbind", data_object$y_list))
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
                                 dist_matrix_list, randcov_params_val) {

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
  w_and_H <- get_w_and_H_spgautor(data_object, dispersion,
                               SigInv, SigInv_X, cov_betahat, Xt_SigInv_X, estmethod)

  w <- w_and_H$w
  # H <- w_and_H$H
  mHldet <- w_and_H$mHldet

  betahat <- tcrossprod(cov_betahat, SigInv_X) %*% w
  X <- data_object$X
  r <- w - X %*% betahat
  rt_SigInv_r <- crossprod(r, SigInv) %*% r

  # get wolfinger objects
  y <- data_object$y
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


get_w_and_H_spglm <- function(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod) {

  family <- data_object$family
  SigInv <- Matrix::bdiag(SigInv_list)
  Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
  y <- as.vector(do.call("rbind", data_object$y_list))
  size <- data_object$size
  w <- get_w_init(family, y, dispersion)
  wdiffmax <- Inf
  iter <- 0
  while (iter < 50 && wdiffmax > 1e-4) {


    iter <- iter + 1
    # compute the d vector
    d <-  get_d(family, w, y, size, dispersion)
    # and then the gradient vector
    g <-  d - Ptheta %*% w
    # Next, compute H
    D <- get_D(family, w, y, size, dispersion)
    D_diag <- diag(D)
    D_list <- lapply(split(D_diag, data_object$local_index), function(x) Diagonal(x = x))
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

    # DSigInv_eigen <- lapply(DSigInv_list, function(x) eigen(x)) # not symm PD so must use eigen
    # DSigInv_det <- prod(unlist(lapply(DSigInv_eigen, function(x) prod(x$values))))
    # DSigInv_Inv <- Matrix::bdiag(lapply(DSigInv_chol, function(x) tcrossprod(t(t(x$vectors) * x$values), x$vectors)))
    DSigInv_Inv_list <- lapply(DSigInv_list, function(x) solve(x))
    DSigInv_Inv <- Matrix::bdiag(DSigInv_Inv_list)
    HInv <- smw_HInv(AInv = DSigInv_Inv, U = SigInv_X, CInv = cov_betahat_Inv)
    solveHg <- HInv %*% g
    wnew <- w - solveHg
    # check overshoot on loglik surface
    dnew <- get_d(family, wnew, y, size, dispersion)
    gnew <- dnew - Ptheta %*% wnew
    if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem", call. = FALSE)
    if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
    wdiffmax <- max(abs(wnew - w))
    # update w
    w <- wnew
  }

  mHldet <- smw_mHldet(A_list = DSigInv_list, AInv = DSigInv_Inv, U = SigInv_X, C = cov_betahat, CInv = cov_betahat_Inv)
  list(w = w, H = NULL, mHldet = mHldet)
}

get_w_and_H_spgautor <- function(data_object, dispersion, SigInv, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod) {

  family <- data_object$family
  Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
  y <- data_object$y
  size <- data_object$size
  w <- get_w_init(family, y, dispersion)
  wdiffmax <- Inf
  iter <- 0
  while (iter < 50 && wdiffmax > 1e-4) {


    iter <- iter + 1
    # compute the d vector
    d <-  get_d(family, w, y, size, dispersion)
    # and then the gradient vector
    g <-  d - Ptheta %*% w
    # Next, compute H
    D <- get_D(family, w, y, size, dispersion)
    H <-  D - Ptheta # not PD but -H is
    solveHg <- solve(H, g)
    wnew <- w - solveHg
    # check overshoot on loglik surface
    dnew <- get_d(family, wnew, y, size, dispersion)
    gnew <- dnew - Ptheta %*% wnew
    if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem", call. = FALSE)
    if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
    wdiffmax <- max(abs(wnew - w))
    # update w
    w <- wnew
  }

  mHldet <- as.numeric(determinant(-H, logarithm = TRUE)$modulus)
  list(w = w, H = NULL, mHldet = mHldet)
}

# gradient of w (lowcase d)
get_d <- function(family, w, y, size, dispersion) {

  if (family == "poisson") {
    d <- -exp(w) + y
  } else if (family == "nbinomial") {
    d <- dispersion * (y - exp(w)) / (dispersion + exp(w))
  } else if (family == "binomial") {
    d <- y - size * expit(w)
  } else if (family == "Gamma") {
    # d <- (1 / w - y) * 1 / dispersion
    d <- 1 / dispersion * (y * exp(-w) - 1)
  } else if (family == "inverse.gaussian") {
    d <- 1 / dispersion * (y - exp(w)) / exp(2 * w)
  }
  d
}

# Hessian of w (cap D)
get_D <- function(family, w, y, size, dispersion) {

  w <- as.vector(w)

  if (family == "poisson") {
    D_vec <- -exp(w)
  } else if (family == "nbinomial") {
    D_vec <- - (dispersion * exp(w) * (dispersion + y)) / ((dispersion + exp(w))^2)
  } else if (family == "binomial") {
    D_vec <- - size * expit(w) / (1 + exp(w))
  } else if (family == "Gamma") {
    # D_vec <- -1 / w^2 * 1 / dispersion
    D_vec <- - 1 / dispersion * y * exp(-w)
  } else if (family == "inverse.gaussian") {
    D_vec <- 1 / dispersion * (exp(w) - 2 * y) / exp(2 * w)
  }
  D <- Diagonal(x = D_vec)
}

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
  }
  w_init
}

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
    disp_recip <- 1 / dispersion
    l00 <- -2 * sum(dgamma(y, shape = disp_recip, scale = dispersion * mu, log = TRUE))
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    disp_recip <- 1 / dispersion
    l00 <- -2 * sum((log(disp_recip) - log(2 * pi) - 3 * log(y)) / 2 - (disp_recip * (y - mu)^2 / (2 * y * mu^2)))
  }
  l00
}

smw_HInv <- function(AInv, U, CInv) {
  mid <- CInv + t(U) %*% AInv %*% U
  # if (all(mid == 0)) diag(mid) <- diag(mid) + 1e-4
  AInv - (AInv %*% U) %*% solve(mid) %*% (t(U) %*% AInv)
}

smw_mHldet <- function(A_list, AInv, U, C, CInv) {
  Aldet <- sum(unlist(lapply(A_list, function(x) determinant(x, logarithm = TRUE)$modulus))) # must be positive det for -H
  Cldet <- 2 * sum(log(diag(t(chol(C)))))
  midldet <- determinant(CInv + t(U) %*% AInv %*% U, logarithm = TRUE)$modulus
  as.numeric(Aldet + Cldet + midldet)
}

get_DSigInv <- function(D, SigInv) {
  D - SigInv
}

get_DSigInv_parallel <- function(cluster_list) {
  D <- cluster_list$D
  S <- cluster_list$S
  get_DSigInv(D, S)
}
