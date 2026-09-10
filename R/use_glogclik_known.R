#' Use composite likelihood (Gaussian sv) for estimation with known covariance parameters
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param formula A formula
#' @param data Data
#' @param dist_matrix A distance matrix (Euclidean)
#' @param partitoin_list A partition matrix
#'
#' @return The known covariance parameters
#'
#' @noRd
use_glogclik_known <- function(spcov_initial, data_object, dist_matrix_list, partition_list = NULL) {
  # this mirrors use_glogclik() but skips optim() entirely because the
  # covariance parameters are already fixed (known) rather than estimated --
  # the pairwise distance/residual setup below is only needed to report the
  # composite likelihood loss value at those fixed parameters
  dist_vector_list <- lapply(dist_matrix_list, function(x) {
    x <- as.matrix(x)
    x <- x[upper.tri(x)]
  })
  if (any(unlist(dist_vector_list) == 0)) {
    warning("Zero distances observed between at least one pair. Ignoring pairs. If using splm(), consider a different estimation method.", call. = FALSE)
  }
  # OLS residuals stand in for the (unobservable) mean-zero spatial process;
  # spdist() on the 1-D residual vector reuses the same distance machinery to
  # get every pairwise residual difference (residual_i - residual_j)
  residual_list <- lapply(data_object$obdata_list, function(d) residuals(lm(data_object$formula, data = d)))
  residual_matrix_list <- lapply(residual_list, function(x) spdist(xcoord_val = x))
  residual_vector_list <- lapply(residual_matrix_list, function(x) {
    x <- as.matrix(x)
    x <- x[upper.tri(x)]
  })
  residual_vector <- unlist(residual_vector_list)
  if (!is.null(data_object$partition_list)) {
    # zero out (and later drop) pairs that fall in different partition
    # factor levels, since the composite likelihood only pools pairs within
    # the same partition group
    partition_vector_list <- lapply(data_object$partition_list, function(x) {
      x <- as.matrix(x)
      x <- x[upper.tri(x)]
    })
    dist_vector_list <- mapply(d = dist_vector_list, p = partition_vector_list, function(d, p) d * p, SIMPLIFY = FALSE)
    residual_vector_list <- mapply(r = residual_vector_list, p = partition_vector_list, function(r, p) r * p, SIMPLIFY = FALSE)
  }

  dist_vector <- unlist(dist_vector_list)
  # keep only strictly positive distances -- removes co-located pairs and,
  # when partitioning was applied, the cross-partition pairs zeroed out above
  dist_index <- dist_vector > 0
  dist_vector <- dist_vector[dist_index]
  residual_vector <- unlist(residual_vector_list)
  residual_vector <- residual_vector[dist_index]
  residual_vector2 <- residual_vector^2

  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(class(spcov_initial), spcov_initial$initial)

  # get loss val
  glogclikloss_val <- get_glogclikloss(spcov_params_val, residual_vector2, dist_vector)

  # return parameter values and optim output
  optim_output <- list(
    method = NA, control = NA, value = glogclikloss_val,
    counts = NA, convergence = NA,
    message = NA, hessian = NA
  )
  # returning output
  list(
    spcov_params_val = spcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known)
  )
}
