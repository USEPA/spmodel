#' Find relevant products to use in Gaussian log-likelihood calculations
#'
#' @param spcov_params_val A \code{spcov_params} object
#' @param ... other arguments
#'
#' @return The relevant Gaussian log-likelihood products
#'
#' @noRd
gloglik_products <- function(spcov_params_val, ...) {
  UseMethod("gloglik_products", spcov_params_val)
}
#' @export
gloglik_products.exponential <- function(spcov_params_val, data_object, estmethod,
                                         dist_matrix_list, randcov_params_val, ...) {


  # making a covariance matrix
  cov_matrix_list <- get_cov_matrix_list(spcov_params_val, dist_matrix_list, randcov_params_val,
                                         data_object$randcov_list, data_object$partition_list,
                                         diagtol = data_object$diagtol)


  # cholesky products
  if (data_object$parallel) {
    cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
      cluster_list_element <- list(
        c = cov_matrix_list[[l]],
        x = data_object$X_list[[l]],
        y = data_object$y_list[[l]]
      )
    })
    cholprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_cholprods_parallel)
    names(cholprods_list) <- names(cov_matrix_list)
  } else {
    cholprods_list <- mapply(
      c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
      function(c, x, y) get_cholprods(c, x, y),
      SIMPLIFY = FALSE
    )
  }

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
  ## t(X) %*% sigma_inverse %*% y
  Xt_SigInv_y <- crossprod(SqrtSigInv_X, SqrtSigInv_y)
  ## t(X) %*% sigma_inverse %*% X)^(-1) %*% t(X) %*% sigma_inverse %*% y
  betahat <- cov_betahat %*% Xt_SigInv_y
  ## lower chol %*% (y - X %*% beta) r stands for residual
  SqrtSigInv_r <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
  ## residual %*% sigma_inverse %*% residual
  rt_SigInv_r <- crossprod(SqrtSigInv_r, SqrtSigInv_r)

  # using wolfinger notation
  l1 <- sum(unlist(lapply(cholprods_list, function(x) 2 * sum(log(diag(x$Sig_lowchol))))))
  l2 <- as.numeric(rt_SigInv_r)

  # returning relevant quantities
  if (estmethod == "reml") {
    l3 <- 2 * sum(log(diag(Xt_SigInv_X_upchol)))
    return(list(l1 = l1, l2 = l2, l3 = l3))
  }

  if (estmethod == "ml") {
    return(list(l1 = l1, l2 = l2))
  }
}
#' @export
gloglik_products.spherical <- gloglik_products.exponential
#' @export
gloglik_products.gaussian <- gloglik_products.exponential
#' @export
gloglik_products.triangular <- gloglik_products.exponential
#' @export
gloglik_products.circular <- gloglik_products.exponential
#' @export
gloglik_products.none <- gloglik_products.exponential
#' @export
gloglik_products.cubic <- gloglik_products.exponential
#' @export
gloglik_products.pentaspherical <- gloglik_products.exponential
#' @export
gloglik_products.cosine <- gloglik_products.exponential
#' @export
gloglik_products.wave <- gloglik_products.exponential
#' @export
gloglik_products.jbessel <- gloglik_products.exponential
#' @export
gloglik_products.gravity <- gloglik_products.exponential
#' @export
gloglik_products.rquad <- gloglik_products.exponential
#' @export
gloglik_products.magnetic <- gloglik_products.exponential

#' @export
gloglik_products.matern <- gloglik_products.exponential
#' @export
gloglik_products.cauchy <- gloglik_products.exponential
#' @export
gloglik_products.pexponential <- gloglik_products.exponential

#' @export
gloglik_products.car <- function(spcov_params_val, data_object, estmethod,
                                 dist_matrix_list, randcov_params_val, ...) {
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
  betahat <- cov_betahat %*% crossprod(SigInv_X, data_object$y)
  r <- data_object$y - data_object$X %*% betahat
  rt_SigInv_r <- crossprod(r, SigInv %*% r)

  # using wolfinger notation
  l1 <- Sigldet
  l2 <- rt_SigInv_r

  # returning relevant quantities
  if (estmethod == "reml") {
    l3 <- 2 * sum(log(diag(Xt_SigInv_X_upchol)))
    return(list(l1 = l1, l2 = l2, l3 = l3))
  }

  if (estmethod == "ml") {
    return(list(l1 = l1, l2 = l2))
  }
}
#' @export
gloglik_products.sar <- gloglik_products.car
