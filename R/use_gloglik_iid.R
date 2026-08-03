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
# fast path for the "iid errors, no spatial dependence, no random effects"
# case: the covariance matrix is sigma^2 * I, so REML/ML have closed-form
# solutions via ordinary least squares and no numerical optimization (optim())
# over covariance parameters is needed at all
use_gloglik_iid <- function(spcov_initial, estmethod, data_object, dist_matrix_list) {
  if (inherits(spcov_initial, c("car", "sar"))) {
    X <- data_object$X
    y <- data_object$y
  } else {
    X <- do.call("rbind", data_object$X_list)
    y <- do.call("rbind", data_object$y_list)
  }

  lmod <- lm(data_object$formula, data = data_object$obdata)
  sse <- sum(residuals(lmod)^2)
  Xt_X <- crossprod(X, X)

  # l1, l2, l3 are the three log-determinant/quadratic-form pieces of the
  # -2*loglik formula; with iid errors the correlation matrix is the
  # identity, so its log-determinant (l1) is exactly zero
  l1 <- 0 # sum of the logs of the identity (all ones)
  l2 <- sse
  # l3 <- 2 * sum(log(diag(abs(R))))
  # log|X'X| via the Cholesky factor's diagonal is numerically more stable
  # than computing det(Xt_X) directly
  l3 <- 2 * sum(log(diag(chol(Xt_X))))

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

  # reconcile a genuinely estimated ie with the numerical floor
  # spcov_matrix.*() applies internally when building Sigma -- mirrors the
  # equivalent GLM-side reconciliation, see R/floor_estimated_ie.R
  spcov_params_val <- floor_estimated_ie(spcov_params_val, spcov_initial$is_known, data_object$diagtol)

  # return parameter values and optim output
  optim_output <- known_optim_output_stub(minustwologlik)

  # return list
  list(
    spcov_params_val = spcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known)
  )
}
