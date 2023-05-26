#' Extract fitted model residuals
#'
#' @description Extract residuals from a fitted model object.
#'   \code{resid} is an alias.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param type \code{"response"} for response residuals, \code{"pearson"}
#'   for Pearson residuals, or \code{"standardized"} for standardized residuals.
#'   For [splm()] and [spautor()] fitted model objects, the default is \code{"response"}.
#'   For [spglm()] and [spgautor()] fitted model objects, deviance residuals are
#'   also available (\code{"deviance"}) and are the default residual type.
#' @param ... Other arguments. Not used (needed for generic consistency).
#' @param model A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#'
#' @details The response residuals are taken as the response minus the fitted values
#'   for the response: \eqn{y - X \hat{\beta}}. The Pearson residuals are the
#'   response residuals pre-multiplied by their inverse square root.
#'   The standardized residuals are Pearson residuals divided by the square
#'   root of one minus the leverage (hat) value. The standardized residuals are often used to
#'   check model assumptions, as they have mean zero and variance approximately one.
#'
#'   \code{rstandard()} is an alias for \code{residuals(model, type = "standardized")}.
#'
#' @return The residuals as a numeric vector.
#'
#' @name residuals.spmodel
#' @method residuals splm
#' @order 1
#' @export
#'
#' @seealso [augment.spmodel()] [cooks.distance.spmodel()] [hatvalues.spmodel()] [influence.spmodel()]
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
residuals.splm <- function(object, type = "response", ...) {
  if (type == "response") {
    return(object$residuals$response)
  } else if (type == "pearson") {
    return(object$residuals$pearson)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("residuals must be response or pearson or standardized")
  }
}
#' @rdname residuals.spmodel
#' @method resid splm
#' @order 2
#' @export
resid.splm <- residuals.splm

#' @rdname residuals.spmodel
#' @method residuals spautor
#' @order 4
#' @export
residuals.spautor <- residuals.splm

#' @rdname residuals.spmodel
#' @method resid spautor
#' @order 5
#' @export
resid.spautor <- residuals.spautor

#' @rdname residuals.spmodel
#' @method rstandard splm
#' @order 3
#' @export
rstandard.splm <- function(model, ...) {
  residuals.splm(model, type = "standardized")
}

#' @rdname residuals.spmodel
#' @method rstandard spautor
#' @order 6
#' @export
rstandard.spautor <- rstandard.splm
