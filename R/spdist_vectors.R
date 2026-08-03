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
