#' Extract log-likelihood
#'
#' @description Find the log-likelihood of a fitted model when \code{estmethod}
#'   is \code{"ml"} or \code{"reml"}.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()] where \code{estmethod}
#'   is \code{"ml"} or \code{"reml"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The log-likelihood.
#'
#' @name logLik.spmodel
#' @method logLik splm
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' logLik(spmod)
logLik.splm <- function(object, ...) {
  if (object$estmethod %in% c("reml", "ml")) {
    minus2loglik <- object$optim$value
    loglik <- -1 / 2 * minus2loglik
    return(loglik)
  } else {
    stop("log likelihood is only defined for the reml or ml estimation")
  }
}

#' @rdname logLik.spmodel
#' @method logLik spautor
#' @order 2
#' @export
logLik.spautor <- logLik.splm
