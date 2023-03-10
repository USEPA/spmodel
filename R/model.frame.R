#' Extract the model frame from a fitted model object
#'
#' @description Extract the model frame from a fitted model object.
#'
#' @param formula A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A model frame that contains the variables used by the formula
#'   for the fitted model object.
#'
#' @name model.frame.spmodel
#' @method model.frame splm
#' @export
#'
#' @seealso [stats::model.frame()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' model.frame(spmod)
model.frame.splm <- function(formula, ...) {
  # model.frame(formula(formula, ...), data = formula$data, ...) too much customization
  model.frame(formula(formula), data = formula$obdata, drop.unused.levels = TRUE, na.action = na.omit)
}

#' @rdname model.frame.spmodel
#' @method model.frame spautor
#' @export
model.frame.spautor <- function(formula, ...) {
  model.frame(formula(formula), data = formula$data[formula$observed_index, , drop = FALSE], drop.unused.levels = TRUE, na.action = na.omit)
}
