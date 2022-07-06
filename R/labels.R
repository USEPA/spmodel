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
#' @method labels spmod
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' labels(spmod)
labels.spmod <- function(object, ...) {
  labels(terms(formula(object)))
}
