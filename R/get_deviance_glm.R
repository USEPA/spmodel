#' Compute per-observation deviance residual components for GLM-type models
#'
#' @param family The response family
#' @param y The response vector
#' @param fitted_response Fitted values on the response scale
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#' @param dispersion The dispersion parameter (used for \code{"nbinomial"} and \code{"beta"})
#'
#' @return The deviance contribution of each observation
#'
#' @noRd
get_deviance_glm <- function(family, y, fitted_response, size, dispersion) {
  # each branch computes half the per-observation deviance as the family-specific
  # log-likelihood ratio between the saturated model (fitted = observed) and the
  # fitted model; the final deviance_i doubles this to match the usual definition
  if (family == "poisson") {
    half_deviance_i <- ifelse(y == 0, 0, y * log(y / fitted_response)) - (y - fitted_response)
  } else if (family == "binomial") {
    half_deviance_i <- ifelse(y == 0, 0, y * log(y / fitted_response)) +
      ifelse(size - y == 0, 0, (size - y) * log((size - y) / (size - fitted_response)))
  } else if (family == "nbinomial") {
    # hand derived
    half_deviance_i <- ifelse(y == 0, 0, y * (log(y / (y + dispersion)) - log(fitted_response / (fitted_response + dispersion)))) +
      dispersion * (log(fitted_response + dispersion) - log(y + dispersion))
  } else if (family == "Gamma") {
    half_deviance_i <- -log(y / fitted_response) + (y - fitted_response) / fitted_response
  } else if (family == "inverse.gaussian") {
    half_deviance_i <- 0.5 * (y - fitted_response)^2 / (y * fitted_response^2)
  } else if (family == "beta") {
    # lgamma() is used instead of log(gamma()) to avoid NA/overflow for large dispersion
    # constant collects the beta-density normalizing terms (log of gamma function ratios)
    constant <- lgamma(fitted_response * dispersion) + lgamma((1 - fitted_response) * dispersion) - lgamma(y * dispersion) - lgamma((1 - y) * dispersion)
    half_deviance_i <- constant + (y - fitted_response) * dispersion * log(y) + ((1 - y) - (1 - fitted_response)) * dispersion * log(1 - y)
  }
  deviance_i <- 2 * half_deviance_i
  as.numeric(deviance_i)
}
