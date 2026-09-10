#' Diagonal distance of a set of locations' bounding box
#'
#' The length of the diagonal of the bounding box, the largest possible
#' pairwise distance between locations. It is more computationally
#' efficient but less precise than finding the maximum distance between
#' two observations in the domain.
#'
#' @param xcoord,ycoord Numeric vectors of x/y coordinates (same length).
#'
#' @return A single number: the bounding box diagonal length. Halve it to
#'   get what the package elsewhere calls \code{max_halfdist}.
#'
#' @noRd
get_bounding_box_dist <- function(xcoord, ycoord) {
  if (length(xcoord) != length(ycoord)) {
    stop("xcoord and ycoord must have the same length.", call. = FALSE)
  }
  x_range <- range(xcoord)
  y_range <- range(ycoord)
  sqrt((x_range[2] - x_range[1])^2 + (y_range[2] - y_range[1])^2)
}
