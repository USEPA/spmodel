get_hatvalues_glm <- function(w, X, data_object, dispersion) {
  # the hat matrix of the whitened residuals
  V <- get_V(w, data_object$family, data_object$size, dispersion)
  SqrtVInv_X <- sqrt(V) * X # same as diag(sqrt(V)) %*% X
  cov_vhat <- chol2inv(chol(crossprod(SqrtVInv_X, SqrtVInv_X)))
  hatvalues <- diag(SqrtVInv_X %*% tcrossprod(cov_vhat, SqrtVInv_X))
  as.numeric(hatvalues)
}

get_V <- function(w, family, size, dispersion) {

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
    V <-  mu / (1 + (mu / dispersion)) # from Ver Hoef and Boveng 2007
  } else if (family == "Gamma") {
    mu <- exp(w)
    V <- mu^2
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    V <- mu^3
  } else if (family == "beta") {
    mu <- exp(w)
    V <- mu * (1 - mu)
  }
  V
}

get_var_y <- function(w, family, size, dispersion) {
  # var(y) = dispersion * var(mu)
  # when dispersion = 1, var(y) = var(mu)
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
