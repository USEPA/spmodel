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
  residuals$pearson^2 * hatvalues / (p * (1 - hatvalues))
}
