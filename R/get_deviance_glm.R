get_deviance_glm <- function(family, y, fitted_link, size, dispersion, offset) {

  if (!is.null(offset)) {
    fitted_link <- fitted_link + offset # undo w = w - offset for deviance to match glm
  }

  fitted_response <- invlink(fitted_link, family, size)

  # faraway p 157
  # y <- pmax(y, 1e-8) # so deviance  Inf is not calculated
  if (family == "poisson") {
    half_deviance_i <- ifelse(y == 0, 0, y * log(y / fitted_response)) - (y - fitted_response)
    # half_deviance_i <- y * pmax(-1e10, log(y / fitted_response)) - (y - fitted_response)
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
    constant <- log(gamma(fitted_response * dispersion)) + log(gamma((1 - fitted_response) * dispersion)) - log(gamma(y * dispersion)) - log(gamma((1 - y) * dispersion))
    half_deviance_i <- constant + (y - fitted_response) * dispersion * log(y) + ((1 - y) - (1 - fitted_response)) * dispersion * log(1 - y)
  }
  deviance_i <- 2 * half_deviance_i
  as.numeric(deviance_i)
}
