#' Make a spatial distance matrix (Euclidean)
#'
#' @param xcoord_val x-coordinate value
#' @param ycoord_val y-coordinate value
#' @param sparse is the matrix sparse?
#'
#' @return A distance matrix (Euclidean)
#'
#' @noRd
spdist <- function(data, xcoord, ycoord, xcoord_val, ycoord_val, sparse = TRUE) {
  if (missing(data)) {
    if (missing(ycoord_val)) {
      spdist_val <- as.matrix(dist(xcoord_val))
    } else {
      spdist_val <- as.matrix(dist(cbind(xcoord_val, ycoord_val)))
    }
  } else {
    if (missing(ycoord)) {
      spdist_val <- as.matrix(dist(data[[xcoord]]))
    } else {
      spdist_val <- as.matrix(dist(cbind(data[[xcoord]], data[[ycoord]])))
    }
  }

  # turn sparse
  if (sparse) {
    spdist_val <- Matrix(spdist_val, sparse = TRUE)
  }
  spdist_val
}
