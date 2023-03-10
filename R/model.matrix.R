#' Extract the model matrix from a fitted model object
#'
#' @description Extract the model matrix (X) from a fitted model object.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The model matrix (of the fixed effects), whose rows represent
#'   observations and whose columns represent explanatory variables corresponding
#'   to each fixed effect.
#'
#' @name model.matrix.spmodel
#' @method model.matrix splm
#' @export
#'
#' @seealso [stats::model.matrix()]
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' model.matrix(spmod)
model.matrix.splm <- function(object, ...) {
  # model.matrix(formula(object, ...), model.frame(object, ...), ...) too much customization
  model.matrix(object$formula, model.frame(object), contrasts = object$contrasts)
}

#' @rdname model.matrix.spmodel
#' @method model.matrix spautor
#' @export
model.matrix.spautor <- model.matrix.splm
