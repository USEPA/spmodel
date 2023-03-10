response_checks_glm <- function(family, y, size) {

  # checks on y
  if (family == "binomial") {

    if (any(size < 1)) {
      stop("All size values must be at least 1.", call. = FALSE)
    }

    if (any(!is.wholenumber(size))) {
      stop("All size values must be a whole number.", call. = FALSE)
    }

    if (any(y < 0)) {
      stop("All response values must be at least 0.", call. = FALSE)
    }

    if (any(!is.wholenumber(y))) {
      stop("All response values must be a whole number.", call. = FALSE)
    }

    if (all(size == 1)) {
      if (!all(y == 0 | y == 1)) {
        stop("All response values must be 0 or 1. 0 indicates a failure and 1 indicates a success.", call. = FALSE)
      }
    }

  } else if (family == "beta") {
    if (any(y <= 0 | y >= 1)) {
      stop("All response values must be greater than 0 and less than 1.", call. = FALSE)
    }
  } else if (family %in% c("poisson", "nbinomial")) {

    if (any(y < 0)) {
      stop("All response values must be at least 0.", call. = FALSE)
    }

    if (any(!is.wholenumber(y))) {
      stop("All response values must be a whole number.", call. = FALSE)
    }

  } else if (family %in% c("Gamma", "inverse.gaussian")) {

    if (any(y <= 0)) {
      stop("All response values must be greater than 0.", call. = FALSE)
    }

  }
}

# check if whole number
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
