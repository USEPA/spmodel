get_coefficients_glm <- function(betahat, spcov_params, dispersion_params, randcov_params) {
  list(fixed = betahat, spcov = spcov_params,
       dispersion = dispersion_params, randcov = randcov_params)
}
