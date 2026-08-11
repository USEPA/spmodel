#' Compute the average pairwise spatial distance among a set of locations
#'
#' The average Euclidean distance across every distinct pair of locations,
#' computed one location ("row") at a time. 
#'
#' @param xcoord,ycoord Numeric vectors of x/y coordinates (same length).
#' @param sample_size If provided, only this many randomly-selected "anchor"
#'   locations have their row-mean computed (each still measured against
#'   all \code{n} locations), rather than every location, tradomg a small
#'   amount of precision for speed when \code{n} itself is very large. The
#'   default (\code{NULL}) uses every location.
#'
#' @return A single number: the average pairwise distance.
#'
#' @noRd
get_avg_dist <- function(xcoord, ycoord, sample_size = NULL) {
  n <- length(xcoord)
  if (length(ycoord) != n) {
    stop("xcoord and ycoord must have the same length.", call. = FALSE)
  }
  if (n < 2) {
    stop("At least two locations are required.", call. = FALSE)
  }

  anchor_idx <- if (is.null(sample_size) || sample_size >= n) {
    seq_len(n)
  } else {
    sample(seq_len(n), sample_size)
  }

  row_means <- unlist(lapply(anchor_idx, function(i) {
    dist_i <- as.numeric(spdist_vectors2(xcoord[i], ycoord[i], xcoord, ycoord, sparse = FALSE))
    mean(dist_i[-i]) # exclude the (always zero) self-distance
  }))

  mean(row_means)
}
