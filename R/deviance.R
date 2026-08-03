#' Fitted model deviance
#'
#' @description Returns the deviance of a fitted model object.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()],
#'   where \code{estmethod} is \code{"ml"} or \code{"reml"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details For objects estimated using \code{"ml"} or \code{"reml"},
#'   the deviance is twice the difference in log-likelihoods between the
#'   saturated (perfect-fit) model and the fitted model.
#'
#' @return The deviance.
#'
#' @name deviance.spmodel
#' @method deviance splm
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' deviance(spmod)
deviance.splm <- function(object, ...) {
  # deviance (twice the saturated-vs-fitted log-likelihood difference) is only
  # meaningful when the model was fit by (RE)ML, since other estimation methods
  # (e.g., sv-wls) don't optimize a likelihood the deviance can be compared against
  if (object$estmethod %in% c("reml", "ml")) {
    deviance <- object$deviance
    return(deviance)
  } else {
    stop("deviance is only defined for reml or ml estimation")
  }
}

#' @rdname deviance.spmodel
#' @method deviance spautor
#' @order 2
#' @export
deviance.spautor <- deviance.splm
