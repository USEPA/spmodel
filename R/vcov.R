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
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' vcov(spmod)
vcov.spmod <- function(object, ...) {
  type <- "fixed"
  if (type == "fixed") {
    return(object$vcov$fixed)
  }
}
