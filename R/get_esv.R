get_esv <- function(residual_vector2, dist_vector, bins, cutoff, formula) {

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

  # set row names to NULL
  # row.names(esv_out) <- NULL

  esv_out
}

get_esv_cloud <- function(residual_vector2, dist_vector, formula) {

  esv_out <- tibble::tibble(dist = dist_vector, gamma = residual_vector2)

  # set row names to NULL
  # row.names(esv_out) <- NULL

  esv_out
}
