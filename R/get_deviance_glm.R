get_deviance_glm <- function(family, y, fitted_val, size, dispersion) {

  # faraway p 157
  # y <- pmax(y, 1e-8) # so deviance  Inf is not calculated
  if (family == "poisson") {
    half_deviance_i <- ifelse(y == 0, 0, y * log(y / fitted_val)) - (y - fitted_val)
    # half_deviance_i <- y * pmax(-1e10, log(y / fitted_val)) - (y - fitted_val)
  } else if (family == "binomial") {
    half_deviance_i <- ifelse(y == 0, 0, y * log(y / fitted_val)) +
      ifelse(size - y == 0, 0, (size - y) * log((size - y) / (size - fitted_val)))
  } else if (family == "nbinomial") {
    # hand derived
    half_deviance_i <- ifelse(y == 0, 0, y * (log(y / (y + dispersion)) - log(fitted_val / (fitted_val + dispersion)))) +
      dispersion * (log(fitted_val + dispersion) - log(y + dispersion))
  } else if (family == "Gamma") {
    half_deviance_i <- -log(y / fitted_val) + (y - fitted_val) / fitted_val
  } else if (family == "inverse.gaussian") {
    half_deviance_i <- 0.5 * (y - fitted_val)^2 / (y * fitted_val^2)
  }
  deviance_i <- 2 * half_deviance_i
  as.numeric(deviance_i)
}
