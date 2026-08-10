#' Compute the cross-distance matrix between two sets of coordinates
#'
#' @param data A data frame with coordinate columns
#' @param data2 A second data frame with coordinate columns
#' @param xcoord The x-coordinate name
#' @param ycoord The y-coordinate name
#' @param dim_coords The number of coordinate dimensions (1 or 2; any other
#'   value is treated as having no meaningful distance, used for the \code{"none"} covariance type)
#' @param sparse Whether to return the result as a sparse \code{Matrix}
#'
#' @return A \code{NROW(data)} by \code{NROW(data2)} matrix of Euclidean
#'   distances between each row of \code{data} and each row of \code{data2}
#'
#' @noRd
spdist_vectors <- function(data, data2, xcoord, ycoord, dim_coords, sparse = TRUE) {
  # a missing coordinate column silently returns NULL from data[[xcoord]]
  # rather than erroring, which outer() then turns into a 0-length dimension
  # instead of a clear failure -- left unchecked, this produces either a
  # dimension mismatch several lines below ("non-conformable arrays") or, if
  # every coordinate column is missing, a degenerate 0-row distance matrix
  # that only fails much later inside an unrelated Cholesky-based linear
  # algebra call ("invalid 'k' argument") -- neither of which points back to
  # the actual problem. Most callers go through get_prediction_object_splm()/
  # _spglm()/predict_block_splm(), which already check newdata via
  # check_newdata_coords() (R/predict_helpers.R) for a friendlier, earlier
  # message; this is the generic backstop for every other caller (fitting,
  # loocv()/kcv()'s local path, conditional(), decorrelate_newdata(), etc.)
  if (dim_coords %in% c(1, 2)) {
    required_coords <- c(xcoord, if (dim_coords == 2) ycoord)
    missing_coords <- unique(c(setdiff(required_coords, names(data)), setdiff(required_coords, names(data2))))
    if (length(missing_coords) > 0) {
      stop(
        "Coordinate column(s) not found: ", paste0("\"", missing_coords, "\"", collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # storing distances
  if (dim_coords == 1) {
    dist_vector <- sqrt(outer(X = data[[xcoord]], Y = data2[[xcoord]], FUN = function(X, Y) (X - Y)^2))
  } else if (dim_coords == 2) { ## finding 2D distance
    dist_vector_x <- outer(X = data[[xcoord]], Y = data2[[xcoord]], FUN = function(X, Y) (X - Y)^2)
    dist_vector_y <- outer(X = data[[ycoord]], Y = data2[[ycoord]], FUN = function(X, Y) (X - Y)^2)
    dist_vector <- sqrt(dist_vector_x + dist_vector_y)
  } else {
    dist_vector <- matrix(Inf, nrow = NROW(data), ncol = NROW(data2)) ## 0D distance (coords not used for "none")
  }
  if (sparse) {
    dist_vector <- Matrix(dist_vector, sparse = TRUE)
  }
  dist_vector
}

spdist_vectors2 <- function(xcoord_val1, ycoord_val1, xcoord_val2, ycoord_val2, dim_coords = 2, sparse = TRUE) {

  # storing distances
  if (dim_coords == 1) {
    dist_vector <- sqrt(outer(X = xcoord_val1, Y = xcoord_val2, FUN = function(X, Y) (X - Y)^2))
  } else if (dim_coords == 2) { ## finding 2D distance
    dist_vector_x <- outer(X = xcoord_val1, Y = xcoord_val2, FUN = function(X, Y) (X - Y)^2)
    dist_vector_y <- outer(X = ycoord_val1, Y = ycoord_val2, FUN = function(X, Y) (X - Y)^2)
    dist_vector <- sqrt(dist_vector_x + dist_vector_y)
  } else {
    dist_vector <- matrix(Inf, nrow = length(xcoord_val1), ncol = length(xcoord_val2)) ## 0D distance (coords not used for "none")
  }
  if (sparse) {
    dist_vector <- Matrix(dist_vector, sparse = TRUE)
  }
  dist_vector
}
