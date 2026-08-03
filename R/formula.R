#' Model formulae
#'
#' Return formula used by a fitted model object.
#'
#' @param x A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
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
# spautor objects store formula the same way splm objects do, so the splm
# method is reused directly instead of writing a near-identical duplicate
formula.spautor <- formula.splm
