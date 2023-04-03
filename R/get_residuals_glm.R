get_residuals_glm <- function(w, y, data_object, deviance_i, hatvalues, dispersion) {

  # add offset
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }

  residuals_response <- y - invlink(w, data_object$family, data_object$size)

  residuals_deviance <- sign(residuals_response) * sqrt(deviance_i)

  residuals_pearson <- residuals_response / sqrt(get_var_y(w, data_object$family, data_object$size, dispersion))

  residuals_standardized <- residuals_deviance / sqrt(1 - hatvalues) # (I - H on bottom)
  list(response = as.numeric(residuals_response), deviance = as.numeric(residuals_deviance),
       pearson = as.numeric(residuals_pearson), standardized = as.numeric(residuals_standardized))
}
