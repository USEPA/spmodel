#' Compute Cook's distance
#'
#' @description Compute the Cook's distance for each observation from a fitted
#'   model object.
#'
#' @param model A fitted model object from [splm()] or [spautor()].
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
#' @method cooks.distance splm
#' @export
#'
#' @seealso [hatvalues.splm()] [influence.splm()] [residuals.splm()]
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

#' @rdname cooks.distance.splm
#' @method cooks.distance spautor
#' @export
cooks.distance.spautor <- cooks.distance.splm
