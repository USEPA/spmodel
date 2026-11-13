#' Draw shared random effects without forming an observation covariance matrix
#' @noRd
draw_sprnorm_random <- function(data, randcov, partition_factor, samples) {
  if (!length(randcov)) return(NULL)
  design <- get_randcov_design(data, randcov, partition_factor)
  q <- ncol(design$observed)
  effects <- matrix(rnorm(q * samples), q, samples) * sqrt(design$variance)
  unname(as.matrix(design$observed %*% effects))
}
