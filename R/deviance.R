#' Fitted model deviance
#'
#' @description Returns the deviance of a fitted model object.
#'
#' @param object A fitted model object from [splm()] or [spautor()]
#'   where \code{estmethod} is \code{"ml"} or \code{"reml"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details For \code{spmod} objects estimated using \code{"ml"} or \code{"reml"},
#'   the deviance is \eqn{(y - X \beta)^T V (y - X \beta)} for an inverse
#'   covariance matrix \eqn{V}, analogous
#'   to residual sums of (whitened) squares.
#'
#' @return The deviance.
#'
#' @method deviance spmod
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' deviance(spmod)
deviance.spmod <- function(object, ...) {
  if (object$estmethod %in% c("reml", "ml")) {
    deviance <- object$deviance
    return(deviance)
  } else {
    stop("deviance is only defined for reml or ml estimation")
  }
}
