#' Compute an empirical autocovariance function
#'
#' @param residual_vector2 A vector of squared (or cross) residuals
#' @param dist_vector A vector of pairwise distances corresponding to \code{residual_vector2}
#' @param bins The number of distance bins
#' @param cutoff The maximum distance to consider
#' @param formula A formula (unused; kept for a consistent signature with related functions)
#'
#' @return A \code{tibble} with the average distance, autocovariance, and pair
#'   count for each distance bin
#'
#' @noRd
get_eacf <- function(residual_vector2, dist_vector, bins, cutoff, formula) {
  # bin all pairwise distances into equal-width classes up to the cutoff, then
  # average the (cross-)products within each class -- this is the standard
  # binned empirical estimator, analogous to a semivariogram but for covariance
  # compute semivariogram classes
  dist_classes <- cut(dist_vector, breaks = seq(0, cutoff, length.out = bins + 1))

  # compute squared differences within each class
  acov <- tapply(residual_vector2, dist_classes, function(x) mean(x))

  # compute pairs within each class
  np <- tapply(residual_vector2, dist_classes, length)

  # set as zero if necessary
  np <- ifelse(is.na(np), 0, np)

  # compute average distance within each class
  dist <- tapply(dist_vector, dist_classes, mean)

  # return output
  eacf_out <- tibble::tibble(
    bins = factor(levels(dist_classes), levels = levels(dist_classes)),
    dist = as.numeric(dist),
    acov = as.numeric(acov),
    np = as.numeric(np)
  )

  eacf_out
}

#' Compute an empirical autocovariance cloud (unbinned)
#'
#' @param residual_vector2 A vector of squared (or cross) residuals
#' @param dist_vector A vector of pairwise distances corresponding to \code{residual_vector2}
#'
#' @return A \code{tibble} with one row per pair, giving its distance and autocovariance
#'
#' @noRd
get_eacf_cloud <- function(residual_vector2, dist_vector) {
  # no binning/averaging -- every pair is returned as its own row for plotting
  tibble::tibble(dist = dist_vector, acov = residual_vector2)
}
