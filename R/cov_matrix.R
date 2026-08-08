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
cov_matrix <- function(spcov_params, dist_matrix, randcov_params = NULL, randcov_Zs = NULL, partition_matrix = NULL, M = NULL, diagtol = 0) {
  # spatial: M is only supplied for autoregressive (CAR/SAR) models, which need it
  # to build the neighbor-weighted precision structure; geostatistical models pass
  # M = NULL and build the covariance from distances directly instead
  if (is.null(M)) {
    cov_matrix_val <- spcov_matrix(spcov_params, dist_matrix, diagtol = diagtol)
  } else {
    cov_matrix_val <- spcov_matrix(spcov_params, dist_matrix, M)
  }


  # random effects: additional (independent) variance components are additive on
  # the covariance scale, so they're simply summed in
  if (!is.null(randcov_params)) {
    randcov_matrix_val <- randcov_matrix(randcov_params, randcov_Zs)
    cov_matrix_val <- cov_matrix_val + randcov_matrix_val
  }

  # partitioning: an elementwise (Hadamard) product with a 0/1 partition indicator
  # matrix zeroes out covariance between observations in different partitions,
  # enforcing independence across partitions for local/composite fits
  if (!is.null(partition_matrix)) {
    cov_matrix_val <- cov_matrix_val * partition_matrix
  }

  cov_matrix_val
}


#' Assemble the observed-data covariance matrix from its parts
#'
#' @param dist_matrix A distance matrix, already computed for the rows in \code{obdata}
#'   (and already anisotropy-transformed by the caller, if applicable)
#' @param obdata The observed data rows the covariance matrix is being built for
#' @param spcov_params A \code{spcov_params} object
#' @param randcov_params A \code{randcov_params} object
#' @param random The random effect formula (or \code{NULL})
#' @param partition_factor The partition factor formula (or \code{NULL})
#' @param diagtol A diagonal tolerance value
#' @param xlev_list Optional factor levels to enforce for each random effect term
#'   (used so a covariance matrix built for a data subset, e.g. a big-data local
#'   neighborhood, still spans the same random effect levels as the full data)
#'
#' @return A covariance matrix
#'
#' @details Shared by \code{covmatrix()}'s observed-data covariance and
#'   \code{predict()}'s big-data local-neighborhood branch, which both need
#'   "the covariance matrix for this exact set of rows".
#'
#' @noRd
get_obs_cov_matrix <- function(dist_matrix, obdata, spcov_params, randcov_params = NULL,
                               random = NULL, partition_factor = NULL, diagtol = 0, xlev_list = NULL) {
  if (is.null(random)) {
    randcov_Zs <- NULL
  } else {
    randcov_names <- get_randcov_names(random)
    randcov_Zs <- get_randcov_Zs(obdata, randcov_names, xlev_list = xlev_list)
  }
  partition_matrix_val <- partition_matrix(partition_factor, obdata)
  cov_matrix(spcov_params, dist_matrix, randcov_params, randcov_Zs, partition_matrix_val, diagtol = diagtol)
}

#' Create the cross-covariance matrix between two sets of observations
#'
#' @param spcov_params A \code{spcov_params} object
#' @param dist_matrix_cross A cross-distance matrix between the two sets of observations
#' @param randcov_params A \code{randcov_params} object
#' @param randcov_Zs_cross Cross random effect design matrix products
#' @param partition_matrix_cross Cross partition matrix
#'
#' @return A cross-covariance matrix (with the nugget/independent error excluded,
#'   since two distinct observations never share independent error)
#'
#' @noRd
cov_matrix_cross <- function(spcov_params, dist_matrix_cross, randcov_params = NULL, randcov_Zs_cross = NULL, partition_matrix_cross = NULL) {
  # spatial
  # temporarily zero out the nugget/independent-error variance ("ie") before
  # building the cross-covariance, since it only contributes when two observations
  # are the exact same location -- never true across two distinct data sets
  spcov_params_ie <- spcov_params[["ie"]]
  spcov_params[["ie"]] <- 0
  cov_matrix_cross_val <- spcov_matrix(spcov_params, dist_matrix_cross)

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

cov_matrix2 <- function(spcov_params, dist_matrix, randcov_matrix = NULL, partition_matrix = NULL, M = NULL, diagtol = 0) {

  # spatial
  if (is.null(M)) {
    cov_matrix_val <- spcov_matrix(spcov_params, dist_matrix, diagtol = diagtol)
  } else {
    cov_matrix_val <- spcov_matrix(spcov_params, dist_matrix, M)
  }


  # random effects
  if (!is.null(randcov_matrix)) {
    cov_matrix_val <- cov_matrix_val + randcov_matrix
  }

  # partitioning
  if (!is.null(partition_matrix)) {
    cov_matrix_val <- cov_matrix_val * partition_matrix
  }

  # diag_add <- min(1e-4, 1e-4 * sum(spcov_params[["de"]], randcov_params))
  # diag(cov_matrix_val) <- diag(cov_matrix_val) + diag_add #
  # possibly needed for random effects stability with no ie
  cov_matrix_val
}

