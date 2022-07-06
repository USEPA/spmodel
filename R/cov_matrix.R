#' Create a covariance matrix
#'
#' @param spcov_params A \code{spcov_params} object
#' @param dist_matrix A distance matrix (that has already been transformed for anisotropy)
#' @param randcov_params A \code{randcov_params} object
#' @param randcov_Zs Random effect design matrices
#' @param partition_matrix Partition matrix
#' @param M An M matrix for autoregressive models
#'
#' @return A covariance matrix
#'
#' @noRd
cov_matrix <- function(spcov_params, dist_matrix, randcov_params = NULL, randcov_Zs = NULL, partition_matrix = NULL, M = NULL) {

  # spatial
  if (is.null(M)) {
    cov_matrix_val <- spcov_matrix(spcov_params, dist_matrix)
  } else {
    cov_matrix_val <- spcov_matrix(spcov_params, dist_matrix, M)
  }


  # random effects
  if (!is.null(randcov_params)) {
    randcov_matrix_val <- randcov_matrix(randcov_params, randcov_Zs)
    cov_matrix_val <- cov_matrix_val + randcov_matrix_val
  }

  # partitioning
  if (!is.null(partition_matrix)) {
    cov_matrix_val <- cov_matrix_val * partition_matrix
  }

  # diag_add <- min(1e-4, 1e-4 * sum(spcov_params[["de"]], randcov_params))
  # diag(cov_matrix_val) <- diag(cov_matrix_val) + diag_add # needed for random effects stability with no ie
  cov_matrix_val
}


cov_matrix_cross <- function(spcov_params, dist_matrix_cross, randcov_params = NULL, randcov_Zs_cross = NULL, partition_matrix_cross = NULL) {

  # spatial
  spcov_params_ie <- spcov_params[["ie"]]
  # this step is necessary because we use diag() in spcov_matrix, which does not translate
  # well to cov(a, b) as it was built for efficiency with cov(a, a)
  spcov_params[["ie"]] <- 0
  cov_matrix_cross_val <- spcov_matrix(spcov_params, dist_matrix_cross)
  # this step is required because ifelse for Matrix objects coerces to a vector
  # but for matrix objects it remains a matrix. So we have to turn it into a matrix
  # object, do ifelse, and then put it back together as a Matrix object

  ## edit 3/31 this ifelse is no longer required as I realized that only the diagonal
  ## should have ie values. Even if two observations share the same location, they will
  ## have different ie realizations. Because this is a cross-covariance, ie values
  ## are not relevant
  # distzero_ie_matrix_val <- Matrix(ifelse(as.matrix(dist_matrix_cross) == 0, spcov_params_ie, 0), sparse = TRUE)
  # cov_matrix_cross_val <-  cov_matrix_cross_val + distzero_ie_matrix_val

  # random effects
  if (!is.null(randcov_params)) {
    randcov_matrix_cross_val <- randcov_matrix(randcov_params, randcov_Zs_cross)
    cov_matrix_cross_val <- cov_matrix_cross_val + randcov_matrix_cross_val
  }

  # partitioning
  if (!is.null(partition_matrix_cross)) {
    cov_matrix_cross_val <- cov_matrix_cross_val * partition_matrix_cross
  }

  cov_matrix_cross_val
}
