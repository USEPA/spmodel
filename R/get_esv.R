#' Compute an empirical semivariogram
#'
#' @param residual_vector2 A vector of squared residual differences
#' @param dist_vector A vector of pairwise distances corresponding to \code{residual_vector2}
#' @param bins The number of distance bins
#' @param cutoff The maximum distance to consider
#' @param formula A formula (unused; kept for a consistent signature with related functions)
#'
#' @return A \code{tibble} with the average distance, semivariance, and pair
#'   count for each distance bin
#'
#' @noRd
get_esv <- function(residual_vector2, dist_vector, bins, cutoff, formula) {
  # classical (Matheron) empirical semivariogram: bin pairwise distances up to
  # the cutoff, then average squared differences within each bin and halve --
  # gamma(h) = 0.5 * mean((z_i - z_j)^2) for pairs whose distance falls in bin h
  # compute semivariogram classes
  dist_classes <- cut(dist_vector, breaks = seq(0, cutoff, length.out = bins + 1))

  # compute squared differences within each class
  gamma <- tapply(residual_vector2, dist_classes, function(x) mean(x) / 2)

  # compute pairs within each class
  np <- tapply(residual_vector2, dist_classes, length)

  # set as zero if necessary
  np <- ifelse(is.na(np), 0, np)

  # compute average distance within each class
  dist <- tapply(dist_vector, dist_classes, mean)

  # return output
  esv_out <- tibble::tibble(
    bins = factor(levels(dist_classes), levels = levels(dist_classes)),
    dist = as.numeric(dist),
    gamma = as.numeric(gamma),
    np = as.numeric(np)
  )

  esv_out
}

#' Compute Cressie's robust empirical semivariogram
#'
#' @param residual_vector12 A vector of absolute-square-root residual differences
#' @param dist_vector A vector of pairwise distances corresponding to \code{residual_vector12}
#' @param bins The number of distance bins
#' @param cutoff The maximum distance to consider
#' @param formula A formula (unused; kept for a consistent signature with related functions)
#'
#' @return A \code{tibble} with the average distance, robust semivariance, and
#'   pair count for each distance bin
#'
#' @noRd
get_esv_robust <- function(residual_vector12, dist_vector, bins, cutoff, formula) {
  # compute semivariogram classes
  dist_classes <- cut(dist_vector, breaks = seq(0, cutoff, length.out = bins + 1))

  # Cressie's robust estimator: averages sqrt(|differences|) instead of squared
  # differences (residual_vector12 is already on that scale), then raises back
  # to the 4th power with a bias correction -- less sensitive to outlier pairs
  # than the classical estimator above
  # compute squared differences within each class
  gamma <- tapply(residual_vector12, dist_classes, function(x) {
    1 / (0.914 + (0.988 / length(x))) * (mean(x)^4)
  })

  # compute pairs within each class
  np <- tapply(residual_vector12, dist_classes, length)

  # set as zero if necessary
  np <- ifelse(is.na(np), 0, np)

  # compute average distance within each class
  dist <- tapply(dist_vector, dist_classes, mean)

  # return output
  esv_out <- tibble::tibble(
    bins = factor(levels(dist_classes), levels = levels(dist_classes)),
    dist = as.numeric(dist),
    gamma = as.numeric(gamma),
    np = as.numeric(np)
  )

  esv_out
}

#' Compute an empirical semivariogram cloud (unbinned)
#'
#' @param residual_vector2 A vector of squared residual differences
#' @param dist_vector A vector of pairwise distances corresponding to \code{residual_vector2}
#'
#' @return A \code{tibble} with one row per pair, giving its distance and semivariance
#'
#' @noRd
get_esv_cloud <- function(residual_vector2, dist_vector) {
  # no binning/averaging -- every pair is returned as its own row for plotting
  tibble::tibble(dist = dist_vector, gamma = residual_vector2 / 2)
}
