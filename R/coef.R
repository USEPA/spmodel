#' Extract fitted model coefficients
#'
#' @description coef extracts fitted model coefficients from [splm()] or [spautor()]
#'   fitted model objects. \code{coefficients} is an alias for it.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param type \code{"fixed"} for fixed effect coefficients, \code{"spcov"} for
#'   spatial covariance parameter coefficients, or \code{"randcov"} for random effect
#'   variance coefficients. Defaults to \code{"fixed"}. If \code{type = "spcov"}, the
#'   coefficient vector is an [spcov_params()] object (which means that has class
#'   matching the spatial covariance function used).
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A named vector of coefficients.
#'
#' @method coef spmod
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' coef(spmod)
#' coefficients(spmod)
#' coef(spmod, type = "spcov")
coef.spmod <- function(object, type = "fixed", ...) {
  if (type == "fixed") {
    return(object$coefficients$fixed)
  } else if (type == "spcov") {
    spcov_coef <- object$coefficients$spcov
    return(spcov_coef)
  } else if (type == "randcov") {
    return(object$coefficients$randcov)
  } else {
    stop("Invalid type argument. The type argument must be \"fixed\", \"spcov\", or \"randcov\".", call. = FALSE)
  }
}
#' @rdname coef.spmod
#' @export
coefficients.spmod <- coef.spmod
