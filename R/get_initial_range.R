#' Get initial range values
#'
#' @param spcov_type The spatial covariance type
#' @param ... max halfdist values (half the maximum distance in the domain)
#'
#' @return The initial range value
#'
#' @noRd
get_initial_range <- function(spcov_type, ...) {
  UseMethod("get_initial_range", structure(list(), class = spcov_type))
}
#' @export
get_initial_range.exponential <- function(spcov_type, max_halfdist, ...) {
  max_halfdist / 3 # effective range of 3
}
#' @export
get_initial_range.spherical <- function(spcov_type, max_halfdist, ...) {
  max_halfdist / 1 # effective range
}
#' @export
get_initial_range.gaussian <- function(spcov_type, max_halfdist, ...) {
  max_halfdist / sqrt(3) # effective range
}
#' @export
get_initial_range.triangular <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.circular <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.none <- function(spcov_type, max_halfdist, ...) {
  Inf
}

#' @export
get_initial_range.ie <- get_initial_range.none

#' @export
get_initial_range.cubic <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.pentaspherical <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.cosine <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.wave <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.jbessel <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.gravity <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.rquad <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.magnetic <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.matern <- function(spcov_type, max_halfdist, ...) {
  max_halfdist / 3 # effective range
}
#' @export
get_initial_range.cauchy <- function(spcov_type, max_halfdist, ...) {
  max_halfdist
}
#' @export
get_initial_range.pexponential <- function(spcov_type, max_halfdist, ...) {
  max_halfdist / 3
}
