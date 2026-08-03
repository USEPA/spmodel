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

  residuals_standardized <- residuals_deviance / sqrt(1 - hatvalues) # (I - H on bottom)
  list(
    response = as.numeric(residuals_response), deviance = as.numeric(residuals_deviance),
    pearson = as.numeric(residuals_pearson), standardized = as.numeric(residuals_standardized)
  )
}
