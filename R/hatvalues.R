#' Compute leverage (hat) values
#'
#' @description Compute the leverage (hat) value for each observation from a fitted
#'   model object.
#'
#' @param model A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Leverage values measure how far an observation's explanatory variables
#'   are relative to the average of the explanatory variables. In other words, observations with high
#'   leverage are typically considered to have an extreme or unusual combination of explanatory
#'   variables. Leverage values are the diagonal of the hat (projection) matrix.
#'   The larger the hat value, the larger the leverage.
#'
#' @return A vector of leverage (hat) values for each observation from the
#'   fitted model object.
#'
#' @method hatvalues splm
#' @export
#'
#' @seealso [cooks.distance.splm()] [influence.splm()] [residuals.splm()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' hatvalues(spmod)
hatvalues.splm <- function(model, ...) {
  model$hatvalues
}

#' @rdname hatvalues.splm
#' @method hatvalues spautor
#' @export
hatvalues.spautor <- hatvalues.splm
