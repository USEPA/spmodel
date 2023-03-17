get_cov_matrix_list <- function(spcov_params, dist_matrix_list, randcov_params, randcov_list, partition_list, diagtol = 0) {
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
