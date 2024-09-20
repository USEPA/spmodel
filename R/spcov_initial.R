#' Create a spatial covariance parameter initial object
#'
#' @description Create a spatial covariance parameter initial object that specifies
#'   initial and/or known values to use while estimating spatial covariance parameters
#'   with [splm()], [spglm()], [spautor()], or [spgautor()].
#'
#' @param spcov_type The spatial covariance function type. Available options include
#'   \code{"exponential"}, \code{"spherical"}, \code{"gaussian"},
#'   \code{"triangular"}, \code{"circular"}, \code{"cubic"},
#'   \code{"pentaspherical"}, \code{"cosine"}, \code{"wave"},
#'   \code{"jbessel"}, \code{"gravity"}, \code{"rquad"},
#'   \code{"magnetic"}, \code{"matern"}, \code{"cauchy"}, \code{"pexponential"},
#'   \code{"car"}, \code{"sar"}, and \code{"none"}.
#' @param de The spatially dependent (correlated) random error variance. Commonly referred to as
#'   a partial sill.
#' @param ie The spatially independent (uncorrelated) random error variance. Commonly referred to as
#'   a nugget.
#' @param range The correlation parameter.
#' @param extra An extra covariance parameter used when \code{spcov_type} is
#'   \code{"matern"}, \code{"cauchy"}, \code{"pexponential"}, \code{"car"}, or
#'   \code{"sar"}.
#' @param rotate Anisotropy rotation parameter (from 0 to \eqn{\pi} radians).
#'   Not used if \code{spcov_type} is \code{"car"} or \code{"sar"}.
#' @param scale Anisotropy scale parameter (from 0 to 1).
#'   Not used if \code{spcov_type} is \code{"car"} or \code{"sar"}.
#' @param known A character vector indicating which spatial covariance parameters are to be
#'   assumed known. The value \code{"given"} is shorthand for assuming all
#'   spatial covariance parameters given to \code{spcov_initial()} are assumed known.
#'
#' @details The \code{spcov_initial} list is later passed to [splm()], [spglm()], [spautor()], or [spgautor()].
#'   \code{NA} values can be given for \code{ie}, \code{rotate}, and \code{scale}, which lets
#'   these functions find initial values for parameters that are sometimes
#'   otherwise assumed known (e.g., \code{rotate} and \code{scale} with [splm()] and [spglm()]
#'   and \code{ie} with [spautor()] and [spgautor()]).
#'   The spatial covariance functions can be generally expressed as
#'   \eqn{de * R + ie * I}, where \eqn{de} is \code{de} above, \eqn{R}
#'   is a matrix that controls the spatial dependence structure among observations,
#'   \eqn{h}, \eqn{ie} is \code{ie} above, and \eqn{I} is and identity matrix.
#'   Note that \eqn{de} and \eqn{ie} must be non-negative while \eqn{range}
#'   must be positive, except when \code{spcov_type} is \code{car} or \code{sar},
#'   in which case \eqn{range} must be between the reciprocal of the maximum
#'   eigenvalue of \code{W} and the reciprocal of the minimum eigenvalue of
#'   \code{W}. Parametric forms for \eqn{R} are given below, where \eqn{\eta = h / range}:
#'   \itemize{
#'     \item exponential: \eqn{exp(- \eta )}
#'     \item spherical: \eqn{(1 - 1.5\eta + 0.5\eta^3) * I(h <= range)}
#'     \item gaussian: \eqn{exp(- \eta^2 )}
#'     \item triangular: \eqn{(1 - \eta) * I(h <= range)}
#'     \item circular: \eqn{(1 - (2 / \pi) * (m * sqrt(1 - m^2) + sin^{-1}(m))) * I(h <= range), m = min(\eta, 1)}
#'     \item cubic: \eqn{(1 - 7\eta^2 + 8.75\eta^3 - 3.5\eta^5 + 0.75\eta^7) * I(h <= range)}
#'     \item pentaspherical: \eqn{(1 - 1.875\eta + 1.25\eta^3 - 0.375\eta^5) * I(h <= range)}
#'     \item cosine: \eqn{cos(\eta)}
#'     \item wave: \eqn{sin(\eta) / \eta * I(h > 0) + I(h = 0)}
#'     \item jbessel: \eqn{Bj(h * range)}, Bj is Bessel-J function
#'     \item gravity: \eqn{(1 + \eta^2)^{-0.5}}
#'     \item rquad: \eqn{(1 + \eta^2)^{-1}}
#'     \item magnetic: \eqn{(1 + \eta^2)^{-1.5}}
#'     \item matern: \eqn{2^{1 - extra}/ \Gamma(extra) * \alpha^{extra} * Bk(\alpha, extra)}, \eqn{\alpha = (2extra * \eta)^{0.5}}, Bk is Bessel-K function wit  order \eqn{1/5 \le extra \le 5}
#'     \item cauchy: \eqn{(1 + \eta^2)^{-extra}}, \eqn{extra > 0}
#'     \item pexponential: \eqn{exp(h^{extra}/range)}, \eqn{0 < extra \le 2}
#'     \item car: \eqn{(I - range * W)^{-1} * M}, weights matrix \eqn{W},
#'      symmetry condition matrix \eqn{M}, observations with no neighbors
#'      are given a unique variance
#'     parameter called \eqn{extra}, \eqn{extra \ge 0}.
#'     \item sar: \eqn{[(I - range * W)(I - range * W)^T]^{-1}},
#'      weights matrix \eqn{W}, \eqn{^T} indicates matrix transpose,
#'       observations with no neighbors are given a unique variance
#'     parameter called \eqn{extra}, \eqn{extra \ge 0}.
#'     \item none: \eqn{0}
#'   }
#'
#'   All spatial covariance functions are valid in one spatial dimension. All
#'   spatial covariance functions except \code{triangular} and \code{cosine} are
#'   valid in two dimensions. An alias for \code{none} is \code{ie}.
#'
#'   When the spatial covariance function is \code{car} or \code{sar}, \code{extra}
#'   represents the variance parameter for the observations in \code{W} without
#'   at least one neighbor (other than itself) -- these are called unconnected
#'   observations. \code{extra} is only used if there is at least one unconnected
#'   observation.
#'
#' @return A list with two elements: \code{initial} and \code{is_known}.
#'   \code{initial} is a named numeric vector indicating the spatial covariance parameters
#'   with specified initial and/or known values. \code{is_known} is a named
#'   numeric vector indicating whether the spatial covariance parameters in
#'   \code{initial} are known or not. The class of the list
#'   matches the value given to the \code{spcov_type} argument.
#'
#' @export
#'
#' @examples
#' # known de value 1 and initial range value 0.4
#' spcov_initial("exponential", de = 1, range = 0.4, known = c("de"))
#' # known ie value 0 and known range value 1
#' spcov_initial("gaussian", ie = 0, range = 1, known = c("given"))
#' # ie given NA
#' spcov_initial("car", ie = NA)
spcov_initial <- function(spcov_type, de, ie, range, extra, rotate, scale, known) {
  if (missing(spcov_type)) {
    stop("spcov_type must be specified", call. = FALSE)
  } else if (!spcov_type %in% c(
    "exponential", "spherical", "gaussian", "triangular", "circular",
    "none", "ie", "cubic", "pentaspherical", "cosine", "wave", "matern", "car", "sar", "jbessel",
    "gravity", "rquad", "magnetic", "cauchy", "pexponential"
  )) {
    stop(paste(spcov_type, "is not a valid spatial covariance function", sep = " "), call. = FALSE)
  }

  # set defaults
  if (missing(de)) {
    de <- NULL
  }

  if (missing(ie)) {
    ie <- NULL
  }

  if (missing(range)) {
    range <- NULL
  }

  if (missing(extra)) {
    extra <- NULL
  }

  if (missing(rotate)) {
    rotate <- NULL
  }

  if (missing(scale)) {
    scale <- NULL
  }

  # paramter checks
  if (!is.null(de) && !is.na(de) && de < 0) {
    stop("de must be positive", call. = FALSE)
  }
  if (!is.null(ie) && !is.na(ie) && ie < 0) {
    stop("ie must be positive", call. = FALSE)
  }
  if (!is.null(range) && !is.na(range) && range < 0 && !spcov_type %in% c("car", "sar")) {
    stop("range must be positive", call. = FALSE)
  }

  if (spcov_type %in% c("matern", "cauchy", "pexponential") && !is.null(extra) && !is.na(extra) && extra < 0) {
    stop("extra must be positive", call. = FALSE)
  }

  if (!is.null(rotate) && !is.na(rotate) && (rotate < 0 || rotate > pi)) {
    stop("rotate must be in [0, pi].", call. = FALSE)
  }

  if (!is.null(scale) && !is.na(scale) && (scale < 0 || scale > 1)) {
    stop("scale must be in [0, 1].", call. = FALSE)
  }

  if (spcov_type == "matern" && !is.null(extra) && !is.na(extra) && (extra < 1 / 5 || extra > 5)) {
    stop("extra must be between 0.2 and 5", call. = FALSE)
  }

  if (spcov_type == "cauchy" && !is.null(extra) && !is.na(extra) && (extra <= 0)) {
    stop("extra must be positive", call. = FALSE)
  }

  if (spcov_type == "pexponential" && !is.null(extra) && !is.na(extra) && (extra <= 0 || extra > 2)) {
    stop("extra must be positive and no larger than 2", call. = FALSE)
  }

  spcov_params_given <- c(
    de = unname(de),
    ie = unname(ie),
    range = unname(range),
    extra = unname(extra),
    rotate = unname(rotate),
    scale = unname(scale)
  )
  if (missing(known)) {
    is_known <- rep(FALSE, length(spcov_params_given))
  } else {
    if (identical(known, "given")) {
      is_known <- rep(TRUE, length(spcov_params_given))
    } else {
      is_known <- names(spcov_params_given) %in% known
    }
  }
  names(is_known) <- names(spcov_params_given)

  # error if NA and known
  spcov_NA <- which(is.na(spcov_params_given))
  if (any(is_known[spcov_NA])) {
    stop("spcov_initial values cannot be NA and known.", call. = FALSE)
  }

  new_spcov_initial <- structure(list(initial = spcov_params_given, is_known = is_known), class = spcov_type)
  new_spcov_initial
}
