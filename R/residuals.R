#' Extract fitted model residuals
#'
#' @description Extract residuals from a fitted model object.
#'   \code{resid} is an alias.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param type \code{"raw"} for raw residuals, \code{"pearson"}
#'   for Pearson residuals, or \code{"standardized"} for standardized residuals.
#'   The default is \code{"raw"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#' @param model A fitted model object from [splm()] or [spautor()].
#'
#' @details The raw residuals are taken as the response minus the fitted values
#'   for the response: \eqn{y - X \hat{\beta}}. The Pearson residuals are the
#'   raw residuals pre-multiplied by their square (Cholesky) root.
#'   The standardized residuals are Pearson residuals divided by the square
#'   root of one minus the leverage (hat) value. The standardized residuals are often used to
#'   check model assumptions, as they have mean zero and variance approximately one.
#'
#'   \code{rstandard()} is an alias for \code{residuals(model, type = "standardized")}.
#'
#' @return The residuals as a numeric vector.
#'
#' @method residuals spmod
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' residuals(spmod)
#' resid(spmod)
#' residuals(spmod, type = "pearson")
#' residuals(spmod, type = "standardized")
#' rstandard(spmod)
residuals.spmod <- function(object, type = "raw", ...) {
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "pearson") {
    return(object$residuals$pearson)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("residuals must be raw or pearson or standardized")
  }
}
#' @rdname residuals.spmod
#' @export
resid.spmod <- residuals.spmod

#' @rdname residuals.spmod
#' @method rstandard spmod
#' @export
rstandard.spmod <- function(model, ...) {
  residuals.spmod(model, type = "standardized")
}
