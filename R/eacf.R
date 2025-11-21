#' Compute the empirical autocovariance
#'
#' @description Compute the empirical autocovariance (i.e., empirical covariance)
#'   for varying bin sizes and cutoff values.
#'
#' @param formula A formula describing the fixed effect structure.
#' @param data A data frame or \code{sf} object containing the variables in \code{formula}
#'   and geographic information.
#' @param xcoord Name of the variable in \code{data} representing the x-coordinate.
#'   Can be quoted or unquoted. Not required if \code{data} is an \code{sf} object.
#' @param ycoord Name of the variable in \code{data} representing the y-coordinate.
#'   Can be quoted or unquoted. Not required if \code{data} is an \code{sf} object.
#' @param cloud A logical indicating whether the empirical autocovariance should
#'   be summarized by distance class or not. When \code{cloud = FALSE} (the default), pairwise autocovariances
#'   are binned and averaged within distance classes. When \code{cloud} = TRUE,
#'   all pairwise autocovariances and distances are returned (this is known as
#'   the "cloud" autocovariance).
#' @param bins The number of equally spaced bins. The default is 15. Ignored if
#'   \code{cloud = TRUE}.
#' @param cutoff The maximum distance considered.
#'   The default is half the diagonal of the bounding box from the coordinates.
#' @param dist_matrix A distance matrix to be used instead of providing coordinate names.
#' @param partition_factor An optional formula specifying the partition factor.
#'   If specified, autocovariances are only computed for observations sharing the
#'   same level of the partition factor.
#'
#' @details The empirical autocovariance (i.e., empirical covariance) is a tool used to visualize and model
#'   spatial dependence by estimating the semivariance of a process at varying distances.
#'   For a constant-mean process, the
#'   autocovariance at distance \eqn{h} is denoted \eqn{Cov(h)} and defined as
#'   \eqn{Cov(z1, z2)}. Under second-order stationarity,
#'   \eqn{Cov(h) = Cov(0) - \gamma(h)}, where \eqn{gamma(h)} is the semivariance function at distance \code{h}. Typically the residuals from an ordinary
#'   least squares fit defined by \code{formula} are second-order stationary with
#'   mean zero. These residuals are used to compute the empirical autocovariance
#'   At a distance \code{h}, the empirical autocovariance is
#'   \eqn{1/N(h) \sum (r1 \times r2)}, where \eqn{N(h)} is the number of (unique)
#'   pairs in the set of observations whose distance separation is \code{h} and
#'   \code{r1} and \code{r2} are residuals corresponding to observations whose
#'   distance separation is \code{h}. In spmodel, these distance bins actually
#'   contain observations whose distance separation is \code{h +- c},
#'   where \code{c} is a constant determined implicitly by \code{bins}. Typically,
#'   only observations whose distance separation is below some cutoff are used
#'   to compute the empirical semivariogram (this cutoff is determined by \code{cutoff}).
#'
#' @name eacf
#'
#' @return If \code{cloud = FALSE}, a tibble (data.frame) with distance bins
#'   (\code{bins}), the average distance (\code{dist}), the average autocovariance (\code{acov}), and the
#'   number of (unique) pairs (\code{np}). If \code{cloud = TRUE}, a tibble
#'   (data.frame) with distance (\code{dist}) and autocovariance (\code{acov})
#'   for each unique pair.
#'
#' @export
#'
#' @examples
#' eacf(sulfate ~ 1, sulfate)
#' plot(eacf(sulfate ~ 1, sulfate))
eacf <- function(formula, data, xcoord, ycoord, cloud = FALSE, bins = 15, cutoff, dist_matrix, partition_factor) {

  # filter out missing response values
  na_index <- is.na(data[[all.vars(formula)[1]]])
  data <- data[!na_index, , drop = FALSE]
  # finding model frame
  data_model_frame <- model.frame(formula, data, drop.unused.levels = TRUE, na.action = na.pass)
  # model matrix with potential NA
  ob_predictors <- complete.cases(model.matrix(formula, data_model_frame))
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }

  # covert sp to sf
  attr_sp <- attr(class(data), "package")
  if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
    stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
  }

  ## convert sf to data frame (point geometry) (1d objects obsolete)
  ### see if data has sf class
  if (inherits(data, "sf")) {
    data <- suppressWarnings(sf::st_centroid(data))
    data <- sf_to_df(data)
    ### name xcoord ".xcoord" to be used later
    xcoord <- ".xcoord"
    ### name ycoord ".ycoord" to be used later
    ycoord <- ".ycoord"
  }

  # compute spatial distances
  if (missing(dist_matrix)) {
    # non standard evaluation for the x and y coordinates
    xcoord <- substitute(xcoord)
    ycoord <- substitute(ycoord)

    if (missing(xcoord)) {
      stop("The xcoord argument must be specified.", call. = FALSE)
    }

    if (!missing(xcoord)) {
      if (!as.character(xcoord) %in% colnames(data)) {
        stop("The xcoord argument must match the name of a variable in data.", call. = FALSE)
      }
    }

    if (missing(ycoord)) {
      dist_matrix <- spdist(data, xcoord)
    } else {
      if (!as.character(ycoord) %in% colnames(data)) {
        stop("The ycoord argument must match the name of a variable in data.", call. = FALSE)
      }
      dist_matrix <- spdist(data, xcoord, ycoord)
    }
  }


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

  dist_vector <- dist_matrix
  if (any(dist_vector == 0)) {
    dist_index <- dist_vector > 0 & dist_vector <= cutoff
  } else {
    dist_index <- dist_vector <= cutoff
  }

  dist_vector <- dist_vector[dist_index]

  # compute squared differences in the residuals
  lmod <- lm(formula = formula, data = data)
  residuals <- residuals(lmod)
  residual_matrix <- outer(residuals, residuals) # ybar is the fitted mean
  residual_matrix <- residual_matrix[upper.tri(residual_matrix)]
  if (!is.null(partition_factor)) {
    residual_matrix <- residual_matrix * partition_matrix_val
  }

  residual_vector <- residual_matrix
  residual_vector <- residual_vector[dist_index]
  residual_vector2 <- residual_vector # just to align with esv() names


  if (cloud) {
    eacf_out <- get_eacf_cloud(residual_vector2, dist_vector)
  } else {
    eacf_out <- get_eacf(residual_vector2, dist_vector, bins, cutoff)
  }




  # remove NA
  # eacf_out <- na.omit(eacf_out)
  eacf_out <- structure(eacf_out, class = c("eacf", class(eacf_out)), call = match.call(), cloud = cloud)
  eacf_out
}

#' @rdname eacf
#' @method plot eacf
#' @param x An object from \code{eacf()}.
#' @param ... Other arguments passed to other methods.
#' @export
plot.eacf <- function(x, ...) {

  cal <- attr(x, "call")
  if (!is.na(m.f <- match("formula", names(cal)))) {
    cal <- cal[c(1, m.f)]
    names(cal)[2L] <- ""
  }
  cc <- deparse(cal, 80)
  nc <- nchar(cc[1L], "c")
  abbr <- length(cc) > 1 || nc > 75
  sub.caption <- if (abbr) {
    paste(substr(cc[1L], 1L, min(75L, nc)), "...")
  } else {
    cc[1L]
  }

  dotlist <- list(...)
  dotlist <- get_eacf_dotlist_defaults(x, dotlist, cloud = attr(x, "cloud"))

  do.call("plot", c(list(x = x$dist, y = x$acov), dotlist))
  title(sub = sub.caption)
}
