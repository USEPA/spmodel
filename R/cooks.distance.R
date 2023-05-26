#' Compute Cook's distance
#'
#' @description Compute the Cook's distance for each observation from a fitted
#'   model object.
#'
#' @param model A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Cook's distance measures the influence of an observation on a fitted
#'   model object. If an observation is influential, its omission from the data
#'   noticeably impacts parameter estimates. The larger the Cook's distance, the
#'   larger the influence.
#'
#' @return A vector of Cook's distance values for each observation from the
#'   fitted model object.
#'
#' @name cooks.distance.spmodel
#' @method cooks.distance splm
#' @order 1
#' @export
#'
#' @seealso [augment.spmodel()] [hatvalues.spmodel()] [influence.spmodel()] [residuals.spmodel()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' cooks.distance(spmod)
cooks.distance.splm <- function(model, ...) {
  model$cooks_distance
}

#' @rdname cooks.distance.spmodel
#' @method cooks.distance spautor
#' @order 2
#' @export
cooks.distance.spautor <- cooks.distance.splm
