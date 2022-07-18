#' Use Gaussian log-likelihood estimation with iid errors and no random effects
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param estmethod The estimation method (\code{"reml"} or \code{"ml"})
#' @param X Model matrix
#' @param y Response vector
#' @param n Sample size
#' @param p Number of fixed effects
#' @param dist_matrix Distance matrix (Euclidean or neighbors)
#'
#' @return Estimated covariance parameters
#'
#' @noRd
use_gloglik_iid <- function(spcov_initial, estmethod, data_object, dist_matrix_list) {
  if (inherits(spcov_initial, c("car", "sar"))) {
    X <- data_object$X
    y <- data_object$y
  } else {
    X <- do.call("rbind", data_object$X_list)
    y <- do.call("rbind", data_object$y_list)
  }
  # X_qr <- qr(X)
  # R <- qr.R(X_qr)
  # R_Inv <- solve(R)
  # betahat <- R_Inv %*% crossprod(qr.Q(X_qr), y)
  # r <- y - X %*% betahat
  # rt_r <- r^2
  # sse <- sum(rt_r)
  lmod <- lm(data_object$formula, data = data_object$obdata)
  sse <- sum(residuals(lmod)^2)
  RX <- crossprod(X, X)

  # interesting discussion here -- should we return the ml likelihood
  # separate from the reml likelihood? We don't have to, as the
  # beta estimates are the same for all values of sigma2, and using the reml
  # likelihood makes it impossible to compare across fixed effects.
  # But what if we want to compare reml models to an iid reml model? I
  # proceed with separate likelihood here, though I am open to changing this.

  l1 <- 0 # sum of the logs of the identity (all ones)
  l2 <- sse
  # l3 <- 2 * sum(log(diag(abs(R))))
  l3 <- sum(log(diag(abs(RX))))

  if (estmethod == "reml") {
    minustwologlik <- as.numeric(l1 + (data_object$n - data_object$p) * log(l2) + l3 + (data_object$n - data_object$p) * (1 + log(2 * pi / (data_object$n - data_object$p))))
    sigma2 <- sse / (data_object$n - data_object$p)
  } else if (estmethod == "ml") {
    minustwologlik <- as.numeric(l1 + data_object$n * log(l2) + data_object$n * (1 + log(2 * pi / data_object$n)))
    sigma2 <- sse / data_object$n
  }
  spcov_params_val <- spcov_initial$initial
  spcov_params_val[["ie"]] <- sigma2

  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_initial), spcov_orig_val = spcov_params_val)

  # return parameter values and optim output
  optim_output <- list(
    method = NA, control = NA, value = minustwologlik,
    counts = NA, convergence = NA,
    message = NA, hessian = NA
  )

  # return list
  list(
    spcov_params_val = spcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known)
  )
}
