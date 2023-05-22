#' Create a spatial covariance parameter object
#'
#' @description Create a spatial covariance parameter object for use with other
#'   functions.
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
#'   A value of 0 (the default) implies no rotation.
#'   Not used if \code{spcov_type} is \code{"car"} or \code{"sar"}.
#' @param scale Anisotropy scale parameter (from 0 to 1).
#'   A value of 1 (the default) implies no scaling.
#'   Not used if \code{spcov_type} is \code{"car"} or \code{"sar"}.
#'
#' @details
#'   Generally, all arguments to \code{spcov_params} must be specified, though
#'   default arguments are often chosen based on \code{spcov_type}.
#'   When \code{spcov_type} is \code{car} or
#'   \code{sar}, \code{ie} is assumed to be 0 unless specified otherwise.
#'   For full parameterizations of all spatial covariance
#'   functions, see [spcov_initial()].
#'
#' @return A named numeric vector of spatial covariance parameters with class \code{spcov_type}.
#'
#' @export
#'
#' @examples
#' spcov_params("exponential", de = 1, ie = 1, range = 1)
spcov_params <- function(spcov_type, de, ie, range, extra, rotate = 0, scale = 1) {
  if (missing(spcov_type)) {
    stop("spcov_type must be specified", call. = FALSE)
  } else if (!spcov_type %in% c(
    "exponential", "spherical", "gaussian", "triangular", "circular",
    "none", "cubic", "pentaspherical", "cosine", "wave", "matern", "car", "sar", "jbessel",
    "gravity", "rquad", "magnetic", "cauchy", "pexponential"
  )) {
    stop(paste(spcov_type), "is not a valid spatial covariance function.")
  }

  if (missing(spcov_type)) {
    stop("spcov_type must be specified.", call. = FALSE)
  }

  if (spcov_type == "none") {
    de <- 0
    range <- Inf
  }

  if (spcov_type %in% c("car", "sar")) {
    if (missing(extra)) {
      stop("extra must be specified. If there are no unconnected sites in the data, set extra = 0.", call. = FALSE)
    }
    if (missing(ie)) {
      ie <- 0
    }
  }

  # some parameter specification checks
  if (spcov_type != "none" && any(missing(de), missing(ie), missing(range))) {
    stop("de, ie, and range must be specified.", call. = FALSE)
  }

  if (spcov_type == "none" && missing(ie)) {
    stop("ie must be specified.", call. = FALSE)
  }

  if (!(spcov_type %in% c("matern", "cauchy", "pexponential", "car", "sar")) && !missing(extra)) {
    stop("extra cannot be specified for this spatial covariance.", call. = FALSE)
  }
  if (spcov_type %in% c("matern", "cauchy", "pexponential") && missing(extra)) {
    stop("extra must be specified.", call. = FALSE)
  }

  if (!is.na(de) && de < 0) {
    stop("de must be positive.", call. = FALSE)
  }
  if (!is.na(ie) && ie < 0) {
    stop("ie must be positive.", call. = FALSE)
  }
  if (!is.na(range) && range < 0 && !spcov_type %in% c("car", "sar")) {
    stop("range must be positive", call. = FALSE)
  }
  if (spcov_type %in% c("matern", "cauchy", "pexponential") && extra < 0) {
    stop("extra must be positive", call. = FALSE)
  }

  if (rotate < 0 || rotate > pi) {
    stop("rotate must be in [0, pi].", call. = FALSE)
  }

  if (scale < 0 || scale > 1) {
    stop("scale must be in [0, 1].", call. = FALSE)
  }

  if (spcov_type == "matern" && (extra < 1 / 5 || extra > 5)) {
    stop("extra must be between 0.2 and 5.", call. = FALSE)
  }

  if (spcov_type == "cauchy" && (extra <= 0)) {
    stop("extra must be positive.", call. = FALSE)
  }

  if (spcov_type == "pexponential" && (extra <= 0 || extra > 2)) {
    stop("extra must be positive and no larger than 2.", call. = FALSE)
  }

  if (spcov_type %in% c("exponential", "spherical", "gaussian", "triangular", "circular", "none", "cubic", "pentaspherical", "cosine", "wave", "jbessel", "gravity", "rquad", "magnetic")) {
    extra <- NULL
  }
  if (spcov_type %in% c("car", "sar")) {
    rotate <- NULL
  }
  if (spcov_type %in% c("car", "sar")) {
    scale <- NULL
  }

  spcov_params_val <- c(
    de = unname(de), ie = unname(ie), range = unname(range),
    extra = unname(extra), rotate = unname(rotate), scale = unname(scale)
  )

  # the constructor giving the class
  new_spcov_params <- structure(spcov_params_val, class = spcov_type)
  new_spcov_params
}
