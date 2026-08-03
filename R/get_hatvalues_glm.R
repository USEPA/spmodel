#' Compute leverage (hat) values for a fitted GLM-type model
#'
#' @param w The latent (link-scale) predictor vector
#' @param X A model matrix
#' @param data_object The data object
#' @param dispersion The dispersion parameter
#'
#' @return Leverage (hat) values for each observation, capped at 0.999 (with
#'   the excess redistributed proportionally to the remaining hat values) to
#'   guard against numerically unstable standardized residuals
#'
#' @noRd
get_hatvalues_glm <- function(w, X, data_object, dispersion) {
  # the hat matrix of the whitened residuals
  # GLM analog of the Gaussian hat matrix: V plays the role Sigma^{-1} plays in
  # the linear-model case, scaling X by sqrt(variance weight) per observation
  V <- get_V(w, data_object$family, data_object$size, dispersion)
  SqrtVInv_X <- sqrt(V) * X # same as diag(sqrt(V)) %*% X
  # cov(betahat) = (X'VX)^{-1}; only the diagonal (leverage) of the resulting
  # hat matrix is needed, so it's built via tcrossprod/diag rather than
  # materializing the full n x n hat matrix
  cov_vhat <- chol2inv(chol(Matrix::forceSymmetric(crossprod(SqrtVInv_X, SqrtVInv_X))))
  hatvalues <- diag(SqrtVInv_X %*% tcrossprod(cov_vhat, SqrtVInv_X))
  # cap extreme leverage values (which would make standardized residuals
  # numerically unstable, since they divide by sqrt(1 - hatvalue)); the excess
  # above the cap is redistributed proportionally so the values still sum to
  # roughly the same total leverage
  if (any(hatvalues > 0.999)) {
    hatvalues_sum <- sum(hatvalues)
    hatvalues[hatvalues > 0.999] <- 0.999
    hatvalues <- hatvalues * (hatvalues_sum / sum(hatvalues))
  }
  as.numeric(hatvalues)
}

#' Compute the GLM variance-function weight \code{V} used in the hat matrix
#'
#' @param w The latent (link-scale) predictor vector
#' @param family The response family
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter
#'
#' @return The weight \eqn{V} such that \eqn{cov(\hat\beta) = (X'VX/dispersion)^{-1}}
#'
#' @noRd
get_V <- function(w, family, size, dispersion) {
  # each branch computes mu (mean on the response scale, via the inverse link)
  # and then the corresponding variance-function weight V(mu), family-specific
  # V = (1 / dispersion) * (dmu / deta)^2 * (V(mu))
  # when canonical link used, dmu / deta = dmu / dtheta = V(mu)
  # so V = (1 / dispersion) * V(mu)
  # V is such that cov(betahat) = (XtVX/dispersion)^{-1}
  # and hence V^{-1}dispersion = Sigma used in fitting
  # but V equals var(y) when dispersion is one
  if (family == "poisson") {
    mu <- exp(w)
    V <- mu
  } else if (family == "binomial") {
    mu <- expit(w)
    V <- size * mu * (1 - mu)
  } else if (family == "nbinomial") {
    mu <- exp(w)
    V <- mu / (1 + (mu / dispersion)) # from Ver Hoef and Boveng 2007
  } else if (family == "Gamma") {
    mu <- exp(w)
    V <- mu^2
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    V <- mu^3
  } else if (family == "beta") {
    mu <- expit(w)
    V <- mu * (1 - mu)
  }
  V
}

#' Compute the response-scale variance \code{var(y)} for a GLM-type family
#'
#' @param w The latent (link-scale) predictor vector
#' @param family The response family
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter
#'
#' @return The response-scale variance for each observation
#'
#' @noRd
get_var_y <- function(w, family, size, dispersion) {
  # var(y) = dispersion * var(mu)
  # when dispersion = 1, var(y) = var(mu)
  # poisson/binomial have no separate dispersion scaling (dispersion is fixed
  # at 1 for these families), so var_y is just get_V(); the other families
  # rescale get_V() by a family-specific true-dispersion factor
  if (family == "poisson") {
    var_y <- get_V(w, family, size, dispersion)
  } else if (family == "binomial") {
    var_y <- get_V(w, family, size, dispersion)
  } else if (family == "nbinomial") {
    mu <- exp(w)
    var_y <- mu + mu^2 / dispersion
  } else if (family == "Gamma") {
    dispersion_true <- 1 / dispersion
    var_y <- get_V(w, family, size, dispersion) * dispersion_true
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    dispersion_true <- 1 / (mu * dispersion)
    var_y <- get_V(w, family, size, dispersion) * dispersion_true
  } else if (family == "beta") {
    dispersion_true <- 1 / (1 + dispersion)
    var_y <- get_V(w, family, size, dispersion) * dispersion_true
  }
  var_y
}
