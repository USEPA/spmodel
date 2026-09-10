#' Find the effective (practical) range of a spatial covariance function
#'
#' The "practical range" convention from geostatistics: the distance at
#' which the spatial correlation function drops to (and, for monotone
#' decaying functions, stays below) \code{target} which is 0.05 by default, the
#' standard convention. Compact and monotone functions have closed forms 
#' except for matern, which is numerically solved for. 
#'
#' Three families need special handling because their correlation functions
#' are not monotone decaying:
#' \itemize{
#'   \item \code{"cosine"} never decays at all and it returns to correlation
#'     1 every \code{2 * pi * range}. The value returned is the first
#'     distance at which correlation drops below \code{target}
#'     (\code{range * acos(target)}), which is \strong{not} a true
#'     effective range (correlation is not small at distances beyond it in
#'     general) and a warning is issued.
#'   \item \code{"wave"} oscillates while decaying, with envelope
#'     \code{1 / (dist / range)}. The value returned
#'     (\code{range / target}) is the distance beyond which the envelope
#'     itself guarantees correlation stays under \code{target} as a
#'     conservative bound, not the first crossing (which happens much
#'     sooner, at \code{pi * range}). A warning is issued.
#'   \item \code{"jbessel"}'s distance argument is \code{dist * range}, not
#'     \code{dist / range} the way every other family works and so effective
#'     range is \strong{inversely} related to \code{range} here (larger
#'     range means faster decay). The value returned uses the standard
#'     large-argument asymptotic envelope for the Bessel J0 function,
#'     \code{sqrt(2 / (pi * x))}. A warning is issued.
#' }
#'
#' \code{"none"}/\code{"ie"} (no spatial term) return 0. \code{"car"}/
#' \code{"sar"} (areal/lattice covariances, used only by
#' \code{spautor()}/\code{spgautor()}) have no continuous distance-based
#' range and are not supported so calling this on them is an error.
#'
#' @param spcov_params An \code{\link{spcov_params}} object.
#' @param target The correlation threshold defining the effective range.
#'   The default, 0.05, is the standard "practical range" convention.
#'
#' @return A single number: the effective range, in the same distance units
#'   as the coordinates \code{spcov_params} was fit/specified with.
#'
#' @noRd
get_effective_range <- function(spcov_params, target = 0.05) {
  if (!is.numeric(target) || length(target) != 1 || target <= 0 || target >= 1) {
    stop("target must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  UseMethod("get_effective_range", spcov_params)
}

# shared helper: matern is the only family with no closed-form inverse
# (correlation involves besselK) -- root-find spcov_vector(spcov_params, d) /
# de == target via adaptive-bound uniroot, reusing spcov_vector()'s own
# correlation formula rather than duplicating it
get_effective_range_monotone <- function(spcov_params, target) {
  de <- spcov_params[["de"]]
  range_val <- spcov_params[["range"]]
  f <- function(d) as.numeric(spcov_vector(spcov_params, d)) / de - target
  upper <- 20 * range_val
  while (f(upper) > 0) {
    upper <- 2 * upper
  }
  uniroot(f, lower = 0, upper = upper, tol = .Machine$double.eps^0.5)$root
}

# shared helper: compact support -- correlation is exactly 0 beyond range
get_effective_range_compact <- function(spcov_params, target) {
  spcov_params[["range"]]
}

# shared helper: rho(d) = (1 + (d/range)^2)^(-p) family (gravity/rquad/
# magnetic are the p = 1/2, 1, 3/2 special cases; cauchy leaves p = extra
# free) -- solving (1+(d*/range)^2)^(-p) = target gives d* = range *
# sqrt(target^(-1/p) - 1)
get_effective_range_invpower <- function(spcov_params, target, p) {
  spcov_params[["range"]] * sqrt(target^(-1 / p) - 1)
}

#' @export
get_effective_range.exponential <- function(spcov_params, target = 0.05) {
  -log(target) * spcov_params[["range"]]
}
#' @export
get_effective_range.gaussian <- function(spcov_params, target = 0.05) {
  sqrt(-log(target)) * spcov_params[["range"]]
}
#' @export
get_effective_range.gravity <- function(spcov_params, target = 0.05) {
  get_effective_range_invpower(spcov_params, target, p = 1 / 2)
}
#' @export
get_effective_range.rquad <- function(spcov_params, target = 0.05) {
  get_effective_range_invpower(spcov_params, target, p = 1)
}
#' @export
get_effective_range.magnetic <- function(spcov_params, target = 0.05) {
  get_effective_range_invpower(spcov_params, target, p = 3 / 2)
}
#' @export
get_effective_range.cauchy <- function(spcov_params, target = 0.05) {
  get_effective_range_invpower(spcov_params, target, p = spcov_params[["extra"]])
}
#' @export
get_effective_range.pexponential <- function(spcov_params, target = 0.05) {
  (-log(target) * spcov_params[["range"]])^(1 / spcov_params[["extra"]])
}
#' @export
get_effective_range.matern <- function(spcov_params, target = 0.05) get_effective_range_monotone(spcov_params, target)

#' @export
get_effective_range.spherical <- function(spcov_params, target = 0.05) get_effective_range_compact(spcov_params, target)
#' @export
get_effective_range.triangular <- get_effective_range.spherical
#' @export
get_effective_range.circular <- get_effective_range.spherical
#' @export
get_effective_range.cubic <- get_effective_range.spherical
#' @export
get_effective_range.pentaspherical <- get_effective_range.spherical

#' @export
get_effective_range.none <- function(spcov_params, target = 0.05) 0
#' @export
get_effective_range.ie <- get_effective_range.none

#' @export
get_effective_range.car <- function(spcov_params, target = 0.05) {
  stop("get_effective_range() is not defined for \"car\"/\"sar\" (areal) covariance types -- there is no continuous distance-based range.", call. = FALSE)
}
#' @export
get_effective_range.sar <- get_effective_range.car

#' @export
get_effective_range.cosine <- function(spcov_params, target = 0.05) {
  warning("\"cosine\" effective range is not well defined. A potential candidate value is used.", call. = FALSE)
  spcov_params[["range"]] * acos(target)
}

#' @export
get_effective_range.wave <- function(spcov_params, target = 0.05) {
  warning("\"wave\" effective range is not well defined. A potential candidate value is used.", call. = FALSE)
  spcov_params[["range"]] / target
}

#' @export
get_effective_range.jbessel <- function(spcov_params, target = 0.05) {
  warning("\"jbessel\" effective range is not well defined. A potential candidate value is used.", call. = FALSE)
  (2 / (pi * target^2)) / spcov_params[["range"]]
}
