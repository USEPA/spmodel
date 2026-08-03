#' Build a covariance matrix for each partition (big-data local indexing)
#'
#' @param spcov_params A \code{spcov_params} object
#' @param dist_matrix_list A list of distance matrices, one per partition
#' @param randcov_params A \code{randcov_params} object (or \code{NULL})
#' @param randcov_list A list of random effect design matrices, one per partition
#' @param partition_list A list of partition matrices, one per partition
#' @param diagtol A diagonal tolerance value
#'
#' @return A list of covariance matrices, one per partition
#'
#' @noRd
get_cov_matrix_list <- function(spcov_params, dist_matrix_list, randcov_params, randcov_list, partition_list, diagtol = 0) {
  # for big-data "local" fitting, the full covariance matrix is never built;
  # instead each partition gets its own (much smaller) covariance matrix, so
  # branch on which optional pieces (random effects, partition weighting)
  # are present and call cov_matrix() once per partition via mapply
  if (is.null(randcov_params) & is.null(partition_list)) {
    cov_matrix_list <- mapply(d = dist_matrix_list, function(d) cov_matrix(spcov_params, d, diagtol = diagtol), SIMPLIFY = FALSE)
  } else if (!is.null(randcov_params) & is.null(partition_list)) {
    cov_matrix_list <- mapply(
      d = dist_matrix_list, r = randcov_list,
      function(d, r) cov_matrix(spcov_params, d, randcov_params, r, diagtol = diagtol), SIMPLIFY = FALSE
    )
  } else if (is.null(randcov_params) & !is.null(partition_list)) {
    cov_matrix_list <- mapply(
      d = dist_matrix_list, p = partition_list,
      function(d, p) cov_matrix(spcov_params, d, partition_matrix = p, diagtol = diagtol), SIMPLIFY = FALSE
    )
  } else {
    cov_matrix_list <- mapply(
      d = dist_matrix_list, r = randcov_list, p = partition_list,
      function(d, r, p) cov_matrix(spcov_params, d, randcov_params, r, p, diagtol = diagtol),
      SIMPLIFY = FALSE
    )
  }
  cov_matrix_list
}
