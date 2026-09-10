#' Get Cook's distance
#'
#' @param residuals A \code{residuals} object.
#' @param hatvalues Leverage vector
#' @param p Number of fixed effects
#'
#' @return Cook's distance vector
#'
#' @noRd
get_cooks_distance <- function(residuals, hatvalues, p) {
  # standard Cook's distance formula: standardized residual squared, scaled
  # by leverage relative to its complement, normalized by the number of
  # fixed effects -- large values flag observations that strongly influence
  # the fitted coefficients
  residuals$standardized^2 * hatvalues / (p * (1 - hatvalues))
}
