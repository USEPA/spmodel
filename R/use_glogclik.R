#' Use composite likelihood (Gaussian sv) for estimation
#'
#' @param spcov_initial A \code{spcov_initial} object
#' @param formula A formula
#' @param data Data
#' @param dist_matrix A distance matrix (Euclidean)
#' @param partition_matrix A partition matrix
#' @param optim_dotlist An optim dotlist
#'
#' @return The covariance parameter estimates
#'
#' @noRd
use_glogclik <- function(spcov_initial, data_object, dist_matrix_list, partition_list = NULL, optim_dotlist) {
  dist_vector_list <- lapply(dist_matrix_list, function(x) {
    x <- as.matrix(x)
    x <- x[upper.tri(x)]
  })
  if (any(unlist(dist_vector_list) == 0)) {
    warning("Zero distances observed between at least one pair. Ignoring pairs. If using splm(), consider a different estimation method.", call. = FALSE)
  }
  residual_list <- lapply(data_object$obdata_list, function(d) residuals(lm(data_object$formula, data = d)))
  residual_matrix_list <- lapply(residual_list, function(x) spdist(xcoord_val = x))
  residual_vector_list <- lapply(residual_matrix_list, function(x) {
    x <- as.matrix(x)
    x <- x[upper.tri(x)]
  })
  residual_vector <- unlist(residual_vector_list)
  if (!is.null(data_object$partition_list)) {
    partition_vector_list <- lapply(data_object$partition_list, function(x) {
      x <- as.matrix(x)
      x <- x[upper.tri(x)]
    })
    dist_vector_list <- mapply(d = dist_vector_list, p = partition_vector_list, function(d, p) d * p, SIMPLIFY = FALSE)
    residual_vector_list <- mapply(r = residual_vector_list, p = partition_vector_list, function(r, p) r * p, SIMPLIFY = FALSE)
  }

  dist_vector <- unlist(dist_vector_list)
  dist_index <- dist_vector > 0
  dist_vector <- dist_vector[dist_index]
  residual_vector <- unlist(residual_vector_list)
  residual_vector <- residual_vector[dist_index]
  residual_vector2 <- residual_vector^2
  # transforming to optim paramters (log scale)
  spcov_orig2optim_val <- spcov_orig2optim(spcov_initial = spcov_initial, spcov_profiled = FALSE)

  # get optim par
  optim_par <- get_optim_par(spcov_orig2optim_val)

  # check optim dotlist
  optim_dotlist <- check_optim_method(optim_par, optim_dotlist)

  # performing optimization
  optim_output <- do.call("optim", c(
    list(
      par = optim_par,
      fn = glogclik,
      spcov_orig2optim = spcov_orig2optim_val,
      residual_vector2 = residual_vector2,
      dist_vector = dist_vector
    ),
    optim_dotlist
  ))

  # transforming to original scale
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim_val, optim_output$par, spcov_profiled = FALSE)
  # making a covariance parameter vector
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim_val), spcov_orig_val = spcov_orig_val)
  # replace range and extra param

  # return parameter values and optim output
  optim_output <- list(
    method = optim_dotlist$method, control = optim_dotlist$control, value = optim_output$value,
    counts = optim_output$counts, convergence = optim_output$convergence,
    message = optim_output$message, hessian = if (optim_dotlist$hessian) optim_output$hessian else FALSE
  )

  # return list
  list(
    spcov_params_val = spcov_params_val,
    optim_output = optim_output, dist_matrix_list = dist_matrix_list,
    is_known = list(spcov = spcov_initial$is_known)
  )
}
