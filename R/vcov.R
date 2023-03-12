#' Calculate variance-covariance matrix for a fitted model object
#'
#' @description Calculate variance-covariance matrix for a fitted model object.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The variance-covariance matrix of coefficients obtained via \code{coef()}.
#'   Currently, only the variance-covariance matrix of the fixed effects is supported.
#'
#' @name vcov.spmodel
#' @method vcov splm
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' vcov(spmod)
vcov.splm <- function(object, ...) {
  type <- "fixed"
  if (type == "fixed") {
    return(object$vcov$fixed)
  }
}

#' @rdname vcov.spmodel
#' @method vcov spautor
#' @order 2
#' @export
vcov.spautor <- vcov.splm
