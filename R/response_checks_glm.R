#' Validate a response vector against its GLM-type family's support
#'
#' @param family The response family
#' @param y The response vector
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#'
#' @return Nothing; errors if \code{y} (or \code{size}) is outside the
#'   permissible range or is not a whole number where required
#'
#' @noRd
response_checks_glm <- function(family, y, size) {
  # each family's likelihood is only defined on a certain support (e.g. counts
  # must be nonnegative integers); values outside that support would silently
  # produce NaN/-Inf log-likelihoods during optimization, so fail fast here instead
  # checks on y
  if (family == "binomial") {
    # size is the number of trials per observation; must be a positive whole number
    if (any(size < 1)) {
      stop("All size values must be at least 1.", call. = FALSE)
    }

    if (any(!is.wholenumber(size))) {
      stop("All size values must be a whole number.", call. = FALSE)
    }

    # y is the number of successes, bounded between 0 and size (whole number)
    if (any(y < 0)) {
      stop("All response values must be at least 0.", call. = FALSE)
    }

    if (any(!is.wholenumber(y))) {
      stop("All response values must be a whole number.", call. = FALSE)
    }

    # with a single trial per observation, successes/failures must be coded 0/1
    if (all(size == 1)) {
      if (!all(y == 0 | y == 1)) {
        stop("All response values must be 0 or 1. 0 indicates a failure and 1 indicates a success.", call. = FALSE)
      }
    }
  } else if (family == "beta") {
    # beta distribution support is the open interval (0, 1)
    if (any(y <= 0 | y >= 1)) {
      stop("All response values must be greater than 0 and less than 1.", call. = FALSE)
    }
  } else if (family %in% c("poisson", "nbinomial")) {
    # count distributions: nonnegative whole numbers
    if (any(y < 0)) {
      stop("All response values must be at least 0.", call. = FALSE)
    }

    if (any(!is.wholenumber(y))) {
      stop("All response values must be a whole number.", call. = FALSE)
    }
  } else if (family %in% c("Gamma", "inverse.gaussian")) {
    # strictly positive, continuous support
    if (any(y <= 0)) {
      stop("All response values must be greater than 0.", call. = FALSE)
    }
  }
}

#' Check whether values are (within tolerance) whole numbers
#'
#' @param x A numeric vector
#' @param tol The numeric tolerance
#'
#' @return A logical vector, \code{TRUE} where \code{x} is within \code{tol} of an integer
#'
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
