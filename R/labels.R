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
#' @method labels splm
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

#' @rdname labels.splm
#' @method labels spautor
#' @export
labels.spautor <- labels.splm
