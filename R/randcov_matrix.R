#' Create a random effects covariance matrix
#'
#' @param randcov_params A \code{randcov_params} object
#' @param randcov_Zs Random effects design matrices
#'
#' @return A random effects covariance matrix
#'
#' @noRd
randcov_matrix <- function(randcov_params = NULL, randcov_Zs = NULL) {
  if (is.null(randcov_params)) {
    randcov_matrix_val <- NULL
  } else {
    randcov_names <- names(randcov_params)
    # var times ZZt
    randcov_matrices <- lapply(randcov_names, function(x) randcov_params[[x]] * randcov_Zs[[x]][["ZZt"]])
    randcov_matrix_val <- Reduce("+", randcov_matrices)
  }
  randcov_matrix_val
}
