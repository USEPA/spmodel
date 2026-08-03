#' Extract the model matrix from a fitted model object
#'
#' @description Extract the model matrix (X) from a fitted model object.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The model matrix (of the fixed effects), whose rows represent
#'   observations and whose columns represent explanatory variables corresponding
#'   to each fixed effect.
#'
#' @name model.matrix.spmodel
#' @method model.matrix splm
#' @order 1
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
  # passing object$contrasts explicitly (rather than letting model.matrix()
  # pick fresh defaults) ensures factor columns are coded exactly as they
  # were when the model was originally fit
  model.matrix(object$formula, model.frame(object), contrasts = object$contrasts)
}

#' @rdname model.matrix.spmodel
#' @method model.matrix spautor
#' @order 2
#' @export
# spautor's model.frame() method already restricts to observed rows, so the
# same construction logic as splm applies unchanged
model.matrix.spautor <- model.matrix.splm
