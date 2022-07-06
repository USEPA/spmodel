#' Create a covariance vector (used in prediction)
#'
#'
#' @param spcov_params Spatial covariance parameters
#' @param dist_vector Distance vector
#' @param randcov_vector Random effects vector
#' @param partition_vector Partition factor vector
#'
#' @return A covariance vector
#'
#' @noRd
cov_vector <- function(spcov_params, dist_vector, randcov_vector = NULL, partition_vector = NULL) {

  # spatial
  spcov_vector_val <- spcov_vector(spcov_params, dist_vector)

  # random effects
  if (!is.null(randcov_vector)) {
    cov_vector_val <- spcov_vector_val + randcov_vector
  } else {
    cov_vector_val <- spcov_vector_val
  }

  # partitioning
  if (!is.null(partition_vector)) {
    cov_vector_val <- cov_vector_val * partition_vector
  }
  cov_vector_val
}
