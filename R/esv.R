#' Compute the empirical semivariogram
#'
#' @description Compute the empirical semivariogram for varying bin sizes and
#'   cutoff values.
#'
#' @param formula A formula describing the fixed effect structure.
#' @param data A data frame, \code{sf}, or \code{sp} object containing the variables in \code{formula}
#'   and geographic information.
#' @param xcoord Name of the variable in \code{data} representing the x-coordinate.
#'   Can be quoted or unquoted. Not required if \code{data} is an \code{sf} or \code{sp} object.
#' @param ycoord Name of the variable in \code{data} representing the y-coordinate.
#'   Can be quoted or unquoted. Not required if \code{data} is an \code{sf} or \code{sp} object.
#' @param dist_matrix A distance matrix to be used instead of providing coordinate names.
#' @param bins The number of equally spaced bins. The default is 15.
#' @param cutoff The maximum distance considered.
#'   The default is half the diagonal of the bounding box from the coordinates.
#' @param partition_factor An optional formula specifying the partition factor.
#'   If specified, semivariances are only computed for observations sharing the
#'   same level of the partition factor.
#'
#' @details The empirical semivariogram is a tool used to visualize and model
#'   spatial dependence by estimating the semivariance of a process at varying distances. The
#'   semivariance at distance \eqn{h} is denoted \eqn{\gamma(h)} and defined as
#'   \eqn{0.5 * Var(z1  - z2)}. Under second-order stationarity,
#'   \eqn{\gamma(h) = Cov(0) - Cov(h)}, where \eqn{Cov(h)} is the covariance function at distance \code{h}. Typically the residuals from an ordinary
#'   least squares fit defined by \code{formula} are second-order stationary with
#'   mean zero. These residuals are used to compute the empirical semivariogram.
#'   At a distance \code{h}, the empirical semivariance is
#'   \eqn{1/N(h) \sum (r1 - r2)^2}, where \eqn{N(h)} is the number of (unique)
#'   pairs in the set of observations whose distance separation is \code{h} and
#'   \code{r1} and \code{r2} are residuals corresponding to observations whose
#'   distance separation is \code{h}. In spmodel, these distance bins actually
#'   contain observations whose distance separation is \code{h +- c},
#'   where \code{c} is a constant determined implicitly by \code{bins}. Typically,
#'   only observations whose distance separation is below some cutoff are used
#'   to compute the empirical semivariogram (this cutoff is determined by \code{cutoff}).
#'
#'   When using [splm()] with \code{estmethod} as \code{"sv-wls"}, the empirical
#'   semivariogram is calculated internally and used to estimate spatial
#'   covariance parameters.
#'
#' @return A data frame with distance bins (\code{bins}), the  average distance
#'   (\code{dist}), the semivariance (\code{gamma}), and the
#'   number of (unique) pairs (\code{np}).
#'
#' @export
#'
#' @examples
#' esv(sulfate ~ 1, sulfate)
esv <- function(formula, data, xcoord, ycoord, dist_matrix, bins = 15, cutoff, partition_factor) {

  # covert sp to sf
  attr_sp <- attr(class(data), "package")
  if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
    # if (inherits(data, c("SpatialPointsDataFrame", "SpatialPolygonsDataFrame"))) {
    # if (!requireNamespace("sf", quietly = TRUE)) { # requireNamespace checks if sf is installed
    #   stop("Install the sf R package to use sp objects in splm()", call. = FALSE)
    # } else {
    #   data <- sf::st_as_sf(data)
    # }
    # data <- sf::st_as_sf(data)
    stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
  }

  ## convert sf to data frame (point geometry) (1d objects obsolete)
  ### see if data has sf class
  if (inherits(data, "sf")) {
    data <- sf_to_df(data)
    ### name xcoord "xcoord" to be used later
    xcoord <- "xcoord"
    ### name ycoord "ycoord" to be used later
    ycoord <- "ycoord"
  }

  # compute spatial distances
  if (missing(dist_matrix)) {
    # non standard evaluation for the x and y coordinates
    xcoord <- substitute(xcoord)
    ycoord <- substitute(ycoord)
    if (missing(ycoord)) {
      dist_matrix <- spdist(data, xcoord)
    } else {
      dist_matrix <- spdist(data, xcoord, ycoord)
    }
  }

  # dist_matrix <- triu(dist_matrix, k = 1) # strict upper triangle
  dist_matrix <- as.matrix(dist_matrix)
  dist_matrix <- dist_matrix[upper.tri(dist_matrix)]

  if (any(dist_matrix == 0)) {
    warning("Zero distances observed between at least one pair. Ignoring pairs. If using splm(), consider a different estimation method.", call. = FALSE)
  }

  if (missing(partition_factor)) {
    partition_factor <- NULL
  }

  if (!is.null(partition_factor)) {
    # partition_matrix_val <- triu(partition_matrix(partition_factor, data = data), k = 1)
    partition_matrix_val <- as.matrix(partition_matrix(partition_factor, data = data))
    partition_matrix_val <- partition_matrix_val[upper.tri(partition_matrix_val)]
    dist_matrix <- dist_matrix * partition_matrix_val
  }

  if (missing(cutoff)) {
    cutoff <- NULL
  }
  if (is.null(cutoff)) {
    cutoff <- max(dist_matrix) / 2
  }
  # dist_vector <- dist_matrix@x # store as a vector
  dist_vector <- dist_matrix
  # dist_index <- dist_vector <= cutoff
  if (any(dist_vector == 0)) {
    dist_index <- dist_vector > 0 & dist_vector <= cutoff
  } else {
    dist_index <- dist_vector <= cutoff
  }

  dist_vector <- dist_vector[dist_index]

  # compute squared differences in the residuals
  lmod <- lm(formula = formula, data = data)
  residuals <- residuals(lmod)
  # residual_matrix <- triu(spdist(xcoord_val = residuals), k = 1)
  residual_matrix <- as.matrix(spdist(xcoord_val = residuals))
  residual_matrix <- residual_matrix[upper.tri(residual_matrix)]
  if (!is.null(partition_factor)) {
    residual_matrix <- residual_matrix * partition_matrix_val
  }
  # residual_vector <- residual_matrix@x # store as a vector
  residual_vector <- residual_matrix
  residual_vector <- residual_vector[dist_index]
  residual_vector2 <- residual_vector^2

  # compute semivariogram classes
  dist_classes <- cut(dist_vector, breaks = seq(0, cutoff, length.out = bins + 1))

  # browser()
  # compute squared differences within each class
  gamma <- tapply(residual_vector2, dist_classes, function(x) mean(x) / 2)

  # compute pairs within each class
  np <- tapply(residual_vector2, dist_classes, length)

  # set as zero if necessary
  np <- ifelse(is.na(np), 0, np)

  # compute average distance within each class
  dist <- tapply(dist_vector, dist_classes, mean)

  # return output
  esv_out <- data.frame(bins = factor(levels(dist_classes), levels = levels(dist_classes)), dist, gamma, np)

  # set row names to NULL
  row.names(esv_out) <- NULL

  # remove NA
  # esv_out <- na.omit(esv_out)
  return(esv_out)
}
