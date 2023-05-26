get_cooks_distance_glm <- function(residuals, hatvalues, p) {
  residuals$standardized^2 * hatvalues / (p * (1 - hatvalues))
}
