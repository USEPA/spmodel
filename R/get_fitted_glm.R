get_fitted_glm <- function(w, data_object) {


  fitted_link <- as.vector(w)
  names(fitted_link) <- data_object$observed_index
  fitted_response <- invlink(fitted_link, data_object$family, data_object$size)
  names(fitted_response) <- data_object$observed_index
  # call latent link?
  fitted_values <- list(response = fitted_response, link = fitted_link)

}

invlink <- function(fitted_link, family, size) {

  if (family == "poisson") {
    fitted <- exp(fitted_link)
  } else if (family == "binomial") {
    fitted <- size * expit(fitted_link)
  } else if (family == "nbinomial") {
    fitted <- exp(fitted_link)
  } else if (family == "Gamma") {
    # fitted <- 1 / fitted_link
    fitted <- exp(fitted_link)
  } else if (family == "inverse.gaussian") {
    fitted <- exp(fitted_link)
  } else if (family == "beta") {
    fitted <- expit(fitted_link)
  }
  fitted
}
