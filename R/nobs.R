#' Extract the number of observations
#'
#' @description Find the number of observations used for model fitting. This
#'   exists so that generic R functions that call \code{nobs()} (e.g.,
#'   \code{stats::nobs()}) work on \code{spmodel} fitted model objects.
#'   Could alternatively rename object$n as object$nobs and rely on stats::nobs.default.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The number of observations used for modeling (i.e., \code{object$n}).
#'
#' @method nobs splm
#' @noRd
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' nobs(spmod)
nobs.splm <- function(object, ...) {
  object$n
}

#' Extract the number of observations for a fitted \code{spautor()} model object
#'
#' @description Copies \code{object$n} call from \code{nobs.splm()}
#'
#' @param object A fitted model object from [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The number of observations used for modeling (i.e., \code{object$n}).
#'
#' @method nobs spautor
#' @noRd
#' @export
#'
#' @examples
#' sealmod <- spautor(log_trend ~ 1, data = seal, spcov_type = "car")
#' nobs(sealmod)
nobs.spautor <- nobs.splm
