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
