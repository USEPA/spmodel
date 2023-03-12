#' Find labels from object
#'
#' @description Find a suitable set of labels from a fitted model object.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A character vector containing the terms used for the fixed effects
#'   from a fitted model object.
#'
#' @name labels.spmodel
#' @method labels splm
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' labels(spmod)
labels.splm <- function(object, ...) {
  labels(terms(formula(object)))
}

#' @rdname labels.spmodel
#' @method labels spautor
#' @order 2
#' @export
labels.spautor <- labels.splm
