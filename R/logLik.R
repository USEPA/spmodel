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
    # number of estimated parameters
    if (object$estmethod == "ml") {
      n_est_param <- object$npar + object$p
    } else {
      n_est_param <- object$npar
    }
    # lm also returns nall which does help with reml warnings of models
    # with different fixed effect structures
    # still insufficient compared to previous version
    loglik <- structure(loglik, nobs = object$n, df = n_est_param, class = "logLik")
    return(loglik)
  } else {
    stop("logLik is only defined if estmethod is \"ml\" or \"reml\".", call. = FALSE)
  }
}

#' @rdname logLik.spmodel
#' @method logLik spautor
#' @order 2
#' @export
logLik.spautor <- logLik.splm
