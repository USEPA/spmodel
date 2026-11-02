#' Simulate independent random effects and a Vecchia spatial-plus-nugget field
#'
#' Random effects are drawn once per level and partition. Spatial neighbors
#' are selected from earlier sites within the same partition by exhaustive
#' distance or covariance ranking. All simulation columns share the same
#' neighborhoods and conditional covariance calculations.
#'
#' @param object A fitted splm covariance context.
#' @param local_list Resolved simulation order, method, and neighborhood size.
#' @param samples Number of simulations.
#' @return A matrix of centered draws in the original row order.
#' @noRd
get_sprnorm_vecchia <- function(object, local_list, samples) {
  data <- object$obdata
  spcov <- coef(object, type = "spcov")
  ord <- local_list$order
  coords <- as.matrix(data[, c(object$xcoord, object$ycoord), drop = FALSE])
  if (object$anisotropy) {
    transformed <- transform_anis(data, object$xcoord, object$ycoord,
      spcov[["rotate"]], spcov[["scale"]])
    coords <- cbind(transformed$xcoord_val, transformed$ycoord_val)
  }
  coords <- unname(coords[ord, , drop = FALSE])
  partition <- partition_group(object$partition_factor, data)
  if (!is.null(partition)) partition <- partition[ord]
  n <- nrow(data)
  neighbors <- sprnorm_vecchia_neighbors(coords, local_list$size,
    local_list$method, spcov, partition)
  Z <- matrix(rnorm(n * samples), n, samples)
  Y <- matrix(0, n, samples)
  target_var <- spcov_target_var(spcov, object$diagtol)
  for (i in seq_len(n)) {
    pool <- neighbors[[i]]
    if (!length(pool)) {
      Y[i, ] <- sqrt(target_var) * Z[i, ]
      next
    }
    distance <- spdist_vectors2(coords[i, 1], coords[i, 2],
      coords[pool, 1], coords[pool, 2], sparse = FALSE)
    cross_cov <- as.numeric(cov_vector(spcov, distance))
    pool_distance <- spdist_vectors2(coords[pool, 1], coords[pool, 2],
      coords[pool, 1], coords[pool, 2], sparse = FALSE)
    pool_cov <- as.matrix(cov_matrix2(spcov, pool_distance, diagtol = object$diagtol))
    chol_pool <- chol(pool_cov)
    weights <- backsolve(chol_pool, forwardsolve(t(chol_pool), cross_cov))
    variance <- max(target_var - sum(weights * cross_cov), 0)
    center <- as.numeric(crossprod(weights, Y[pool, , drop = FALSE]))
    Y[i, ] <- center + sqrt(variance) * Z[i, ]
  }
  result <- matrix(0, n, samples)
  result[ord, ] <- Y
  random <- draw_sprnorm_random(data, coef(object, type = "randcov"),
    object$partition_factor, samples)
  if (!is.null(random)) result <- result + random
  result
}

#' Find spatial neighbors among preceding sites within each partition
#' @noRd
sprnorm_vecchia_neighbors <- function(coords, size, method, spcov, partition = NULL) {
  n <- nrow(coords)
  neighbors <- rep(list(integer()), n)
  groups <- if (is.null(partition)) list(seq_len(n)) else split(seq_len(n), partition, drop = TRUE)
  covariance_rank <- method == "covariance" && inherits(spcov, c("cosine", "wave", "jbessel"))
  for (group in groups) {
    if (length(group) < 2L) next
    for (j in seq.int(2L, length(group))) {
      i <- group[j]
      pool <- group[seq_len(j - 1L)]
      if (length(pool) > size) {
        distance <- as.numeric(spdist_vectors2(coords[i, 1], coords[i, 2],
          coords[pool, 1], coords[pool, 2], sparse = FALSE))
        score <- if (covariance_rank) -abs(as.numeric(cov_vector(spcov, distance))) else distance
        pool <- pool[order(score, pool)[seq_len(size)]]
      }
      neighbors[[i]] <- pool
    }
  }
  neighbors
}

#' Draw shared random effects without forming an observation covariance matrix
#' @noRd
draw_sprnorm_random <- function(data, randcov, partition_factor, samples) {
  if (!length(randcov)) return(NULL)
  design <- get_randcov_design(data, randcov, partition_factor)
  q <- ncol(design$observed)
  effects <- matrix(rnorm(q * samples), q, samples) * sqrt(design$variance)
  unname(as.matrix(design$observed %*% effects))
}
