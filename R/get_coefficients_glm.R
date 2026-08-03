#' Assemble the coefficients list for a fitted GLM-type model
#'
#' @param betahat Estimated fixed effects
#' @param spcov_params A \code{spcov_params} object
#' @param dispersion_params A \code{dispersion_params} object
#' @param randcov_params A \code{randcov_params} object (or \code{NULL})
#'
#' @return A list with elements \code{fixed}, \code{spcov}, \code{dispersion},
#'   and \code{randcov}, as returned by \code{coef()} with the various
#'   \code{type} options
#'
#' @noRd
get_coefficients_glm <- function(betahat, spcov_params, dispersion_params, randcov_params) {
  # bundle the four coefficient types into one list so coef() can dispatch
  # on its type argument without recomputation; GLM models add a dispersion
  # parameter on top of what get_coefficients() stores for splm/spautor
  list(
    fixed = betahat, spcov = spcov_params,
    dispersion = dispersion_params, randcov = randcov_params
  )
}
