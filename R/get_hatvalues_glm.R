#' Compute leverage (hat) values for a fitted GLM-type model
#'
#' @param w The latent (link-scale) predictor vector (without the offset added)
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
  # w excludes the offset up to this point,
  # so add it back in before evaluating the leverage weight. The weight below
  # is a function of the fitted mean, and the offset is part of the linear
  # predictor that produces that mean. An observation with a large offset has
  # a large fitted mean independent of the covariates. Leaving the
  # offset out evaluates the weight at the wrong mean, which tilts leverage
  # toward the wrong observations; stats::glm() likewise builds its hat matrix
  # from working weights evaluated at the offset-inclusive linear predictor.
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }

  # the hat matrix of the whitened residuals
  # GLM analog of the Gaussian hat matrix: V plays the role Sigma^{-1} plays in
  # the linear-model case, scaling X by sqrt(variance weight) per observation
  V <- get_V(w, data_object$family, data_object$size, dispersion)
  SqrtVInv_X <- sqrt(V) * X # same as diag(sqrt(V)) %*% X
  # cov(betahat) = (X'VX)^{-1}; only the diagonal (leverage) of the resulting
  # hat matrix is needed, so get_diag_XVXt() forms it without ever
  # materializing the full n x n hat matrix
  cov_vhat <- chol2inv(chol(Matrix::forceSymmetric(crossprod(SqrtVInv_X, SqrtVInv_X))))
  hatvalues <- get_diag_XVXt(SqrtVInv_X, cov_vhat)
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

#' Compute the leverage weight \code{V} used in the GLM hat matrix
#'
#' @param w The latent (link-scale) predictor vector, with any offset already
#'   added
#' @param family The response family
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter
#'
#' @return The per-observation expected information on the link (\code{w})
#'   scale, \eqn{E(-D_{ii})}, for use as the weight in the hat matrix
#'
#' @noRd
get_V <- function(w, family, size, dispersion) {
  # Leverage is a projection carried out for whichever metric the estimator
  # actually corresponds to, so the weight here is the expected information each
  # observation contributes on the w (link) scale: E(-D_ii), i.e. the
  # expectation over y of the D_ii expressions in get_D(). Equivalently, this
  # is the iteratively reweighted least squares working weight
  # (dmu / deta)^2 / V(mu) that stats::glm() uses to build its hat matrix.
  #
  # This is NOT in general the variance function V(mu), and it is never var(y).
  # For an exponential dispersion family the identity is
  #
  #   E(-D_ii) = (dmu / deta)^2 / (a(phi) * V(mu))
  #
  # which collapses to V(mu) / a(phi) only under a canonical link. spmodel uses
  # log and logit links throughout, which are canonical for the Poisson and the
  # binomial but not for the gamma, inverse Gaussian, negative binomial, or
  # beta. The negative binomial branch below was already the correct weight and
  # is unchanged; the gamma and inverse Gaussian previously used V(mu), which
  # weighted observations by a factor varying with mu when it should not.
  #
  # Any constant shared by all observations cancels out of the hat matrix
  # so the constant weights returned below for the
  # gamma and inverse Gaussian give the ordinary least squares hat matrix. That
  # is the same answer glm(family = Gamma(link = "log")) returns.
  w <- as.vector(w)
  n <- length(w)
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
    # D_ii = -dispersion * y / mu, so E(-D_ii) = dispersion * E(y) / mu =
    # dispersion, free of mu
    V <- rep(dispersion, n)
  } else if (family == "inverse.gaussian") {
    # D_ii = -dispersion * (mu^2 + y^2) / (2 * y * mu). With shape
    # lambda = mu * dispersion, E(1 / y) = 1 / mu + 1 / (mu * dispersion), so
    # E(-D_ii) = dispersion + 1 / 2, again free of mu
    V <- rep(dispersion + 0.5, n)
  } else if (family == "beta") {
    mu <- expit(w)
    # the -2 * sinh(w) * k0 half of k1 in get_D() has expectation zero, since
    # E(k0) = 0, leaving the trigamma term
    V <- dispersion^2 * (mu * (1 - mu))^2 *
      (trigamma(mu * dispersion) + trigamma((1 - mu) * dispersion))
  }
  V
}

#' Compute the response-scale variance \code{var(y)} for a GLM-type family
#'
#' @param w The latent (link-scale) predictor vector, with any offset already
#'   added
#' @param family The response family
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter
#'
#' @return The response-scale variance \eqn{var(y_i)} for each observation
#'
#' @noRd
get_var_y <- function(w, family, size, dispersion) {
  # The actual variance of the response, used to standardize Pearson residuals
  # so that they have variance approximately one.
  #
  # This is a different quantity from get_V(), and the two are no longer related
  # by a rescaling, so each is computed directly here rather than one in terms
  # of the other. get_V() returns the link-scale information weight used for
  # leverage, which carries no dispersion (any constant cancels out of the hat
  # matrix); var(y) does carry it.
  #
  # Dividing the Pearson residual by sqrt(var(y)) departs from stats::glm(),
  # which divides by sqrt(V(mu)) alone and deliberately leaves the dispersion in
  # the residual so that sum(r^2) / (n - p) estimates it, then restores it in
  # rstandard.glm(). spmodel estimates the dispersion by (restricted) maximum likelihood
  # instead, so nothing downstream needs the dispersion left in, and residuals
  # on a unit scale are the desired useful diagnostic.
  w <- as.vector(w)
  if (family == "poisson") {
    var_y <- exp(w)
  } else if (family == "binomial") {
    mu <- expit(w)
    var_y <- size * mu * (1 - mu)
  } else if (family == "nbinomial") {
    mu <- exp(w)
    var_y <- mu + mu^2 / dispersion
  } else if (family == "Gamma") {
    var_y <- exp(2 * w) / dispersion # mu^2 / dispersion
  } else if (family == "inverse.gaussian") {
    # shape lambda = mu * dispersion, so var(y) = mu^3 / lambda = mu^2 / dispersion
    var_y <- exp(2 * w) / dispersion
  } else if (family == "beta") {
    mu <- expit(w)
    var_y <- mu * (1 - mu) / (1 + dispersion)
  }
  var_y
}

#' Compute the dispersion factor \eqn{a(\varphi)} for a GLM-type family
#'
#' @param w The latent (link-scale) predictor vector, with any offset already
#'   added
#' @param family The response family
#' @param size Binomial trial sizes (unused; kept for a notational consistency)
#' @param dispersion The dispersion parameter
#'
#' @return The factor \eqn{a_i(\varphi)} relating the deviance contribution
#'   returned by \code{get_deviance_glm()} to the scaled deviance
#'
#' @noRd
get_dispersion_factor <- function(w, family, size, dispersion) {
  # This is the third of three distinct family "variances" used in this file,
  # and they are easy to mix up:
  #
  #   get_V()                  link-scale information weight, for leverage
  #   get_var_y()              var(y_i), for Pearson residuals
  #   get_dispersion_factor()  a_i(phi), for standardized residuals
  #
  # a_i(phi) is defined by var(y) = a_i(phi) * V(mu), and equivalently relates
  # the deviance contribution d_i that get_deviance_glm() returns to the scaled
  # deviance (twice the saturated-minus-fitted log-likelihood ratio):
  #
  #   scaled deviance_i = d_i / a_i(phi)
  #
  # It is important because the scaled deviance is the quantity with expectation
  # approximately one for every family, so dividing the deviance residual by
  # sqrt(a_i * (1 - h_i)) is what makes a standardized residual actually
  # standardized. stats::rstandard.glm() does the same thing.
  #
  # get_deviance_glm() already returns the scaled deviance for the Poisson,
  # binomial, negative binomial, and beta (matching glm(), glm.nb(), and
  # betareg()), so a_i is 1 there and those families are unaffected. The gamma
  # and inverse Gaussian instead use glm()'s unscaled deviance, so they need
  # the factor below.
  #
  # Note that a_i is NOT spmodel's dispersion parameter.
  # spmodel parameterizes the gamma distribution so that
  # var(y) = mu^2 / dispersion, making a_i the reciprocal of dispersion (glm()
  # would report 1 / dispersion). It parameterizes the inverse Gaussian with
  # shape lambda = mu * dispersion, so var(y) = mu^2 / dispersion and a_i
  # varies with the fitted mean rather than being a single number at all.
  w <- as.vector(w)
  n <- length(w)
  if (family %in% c("poisson", "binomial", "nbinomial", "beta")) {
    a_phi <- rep(1, n)
  } else if (family == "Gamma") {
    a_phi <- rep(1 / dispersion, n)
  } else if (family == "inverse.gaussian") {
    a_phi <- 1 / (exp(w) * dispersion) # 1 / (mu * dispersion) = 1 / lambda
  }
  a_phi
}
