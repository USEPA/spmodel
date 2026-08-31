#' Get residuals for a fitted GLM-type model
#'
#' @param w The latent (link-scale) predictor vector (without the offset added)
#' @param y Response vector
#' @param data_object The data object
#' @param deviance_i Per-observation deviance contributions
#' @param hatvalues Leverage values
#' @param dispersion The dispersion parameter
#'
#' @return A list of relevant residuals (response, deviance, Pearson, and standardized)
#'
#' @noRd
get_residuals_glm <- function(w, y, data_object, deviance_i, hatvalues, dispersion) {
  # w excludes the offset up to this point (kept separate for the optimizer),
  # so add it back in before mapping to the response scale
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }

  # invlink() maps the latent (link-scale) predictor back to the response scale
  residuals_response <- y - invlink(w, data_object$family, data_object$size)

  # deviance residuals keep the sign of the raw residual but use the
  # per-observation deviance contribution as their magnitude
  residuals_deviance <- sign(residuals_response) * sqrt(deviance_i)

  # Pearson residuals scale the raw residual by the family's mean-variance relationship
  residuals_pearson <- residuals_response / sqrt(get_var_y(w, data_object$family, data_object$size, dispersion))

  # Standardized (studentized) deviance residuals. 
  #
  # The deviance contribution d_i that get_deviance_glm() returns has
  # expectation approximately a_i(phi).
  # Dividing d_i by a_i first recovers the scaled deviance, whose expectation is
  # approximately one for every family; sqrt(1 - h_i) then corrects for the
  # shrinkage caused by having fitted the mean. Equivalently, divide the deviance residual by sqrt(a_i * (1 - h_i)), which is the
  # same construction stats::rstandard.glm() uses.
  #
  # a_i is 1 for the Poisson, binomial, negative binomial, and beta, so those
  # four families are numerically unchanged; only the gamma and inverse
  # Gaussian are rescaled. Note this makes the standardized residual something
  # other than the deviance residual divided by sqrt(1 - h) for those two
  # families, exactly as it is in rstandard.glm().
  a_phi <- get_dispersion_factor(w, data_object$family, data_object$size, dispersion)
  residuals_standardized <- residuals_deviance / sqrt(a_phi * (1 - hatvalues)) # (I - H on bottom)
  list(
    response = as.numeric(residuals_response), deviance = as.numeric(residuals_deviance),
    pearson = as.numeric(residuals_pearson), standardized = as.numeric(residuals_standardized)
  )
}
