#' Extract the model frame from a fitted model object
#'
#' @description Extract the model frame from a fitted model object.
#'
#' @param formula A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A model frame that contains the variables used by the formula
#'   for the fitted model object.
#'
#' @name model.frame.spmodel
#' @method model.frame splm
#' @order 1
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
  # na.action = na.omit and drop.unused.levels = TRUE are fixed rather than
  # passed through ... so the returned frame always matches the rows/factor
  # levels actually used when the model was fit
  model.frame(formula(formula), data = formula$obdata, drop.unused.levels = TRUE, na.action = na.omit)
}

#' @rdname model.frame.spmodel
#' @method model.frame spautor
#' @order 2
#' @export
model.frame.spautor <- function(formula, ...) {
  # spautor() keeps both observed and missing rows in $data (autoregressive
  # models need the full neighborhood structure), so subset to observed_index
  # here to get just the rows the model was fit to
  model.frame(formula(formula), data = formula$data[formula$observed_index, , drop = FALSE], drop.unused.levels = TRUE, na.action = na.omit)
}
