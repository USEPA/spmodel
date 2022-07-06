#' Extract model fitted values
#'
#' @description Extract fitted values from fitted model objects. \code{fitted.values}
#'   is an alias.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param type \code{"response"} for fitted values of the response, \code{"spcov"}
#'   for fitted values of the spatial random errors, or \code{"randcov"} for
#'   fitted values of the random effects. The default is \code{"response"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details When \code{type} is \code{"response"}, the fitted values
#'   for each observation are the standard fitted values \eqn{X \hat{\beta}}.
#'   When \code{type} is \code{"spcov"} the fitted values for each observation
#'   are (generally) the best linear unbiased predictors of the spatial dependent and spatial
#'   independent random error. When \code{type} is \code{"randcov"}, the fitted
#'   values for each level of each random effect are (generally) the best linear unbiased
#'   predictors of the corresponding random effect. The fitted values for \code{type}
#'   \code{"spcov"} and \code{"randcov"} can be used to assess assumptions
#'   for each component of the fitted model object (e.g., assess a Gaussian assumption).
#'
#' @return The fitted values according to \code{type}.
#'
#' @method fitted spmod
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' fitted(spmod)
#' fitted.values(spmod)
#' fitted(spmod, type = "spcov")
fitted.spmod <- function(object, type = "response", ...) {
  if (type == "response") {
    fitted_val <- object$fitted$response
  } else if (type == "spcov") {
    fitted_val <- object$fitted$spcov
  } else if (type == "randcov") {
    fitted_val <- object$fitted$randcov
  } else {
    stop("Invalid type argument. The type argument must be \"response\", \"spcov\", or \"randcov\".", call. = FALSE)
  }
  fitted_val
}
#' @rdname fitted.spmod
#' @method fitted.values spmod
#' @export
fitted.values.spmod <- fitted.spmod
