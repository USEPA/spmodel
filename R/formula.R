#' Model formulae
#'
#' Return formula used by a fitted model object.
#'
#' @param x A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The formula used by a fitted model object.
#'
#' @method formula spmod
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' formula(spmod)
formula.spmod <- function(x, ...) {
  formula(x$formula)
}
