#' Compute Cook's distance for a fitted GLM-type model
#'
#' @param residuals A \code{data.frame} of residuals with a \code{standardized} column
#' @param hatvalues Leverage (hat) values
#' @param p The number of estimated fixed effects
#'
#' @return Cook's distance for each observation
#'
#' @noRd
get_cooks_distance_glm <- function(residuals, hatvalues, p) {
  # same Cook's distance formula as the linear-model version, applied to
  # GLM standardized residuals and hat values
  residuals$standardized^2 * hatvalues / (p * (1 - hatvalues))
}
