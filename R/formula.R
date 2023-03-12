#' Model formulae
#'
#' Return formula used by a fitted model object.
#'
#' @param x A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The formula used by a fitted model object.
#'
#' @name formula.spmodel
#' @method formula splm
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' formula(spmod)
formula.splm <- function(x, ...) {
  formula(x$formula)
}

#' @rdname formula.spmodel
#' @method formula spautor
#' @order 2
#' @export
formula.spautor <- formula.splm
