#' Regression diagnostics
#'
#' @description Provides basic quantities which are used in forming
#'   a wide variety of diagnostics for checking the quality of fitted model objects.
#'
#' @param model A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details This function calls [residuals.spmodel()], [hatvalues.spmodel()],
#'   and [cooks.distance.spmodel()] and puts the results into a tibble. It is
#'   primarily used when calling [augment.spmodel()].
#'
#' @return A tibble with residuals (\code{.resid}), leverage values (\code{.hat}),
#'   cook's distance (\code{.cooksd}), and standardized residuals (\code{.std.resid}).
#'
#' @name influence.spmodel
#' @method influence splm
#' @order 1
#' @export
#'
#' @seealso [augment.spmodel()] [cooks.distance.spmodel()] [hatvalues.spmodel()] [residuals.spmodel()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' influence(spmod)
influence.splm <- function(model, ...) {
  # standardized residuals are requested via a separate residuals() call
  # (rather than derived here from .resid and .hat) since the standardization
  # accounts for the full model covariance, not just leverage
  tibble::tibble(
    .resid = residuals(model),
    .hat = hatvalues(model),
    .cooksd = cooks.distance(model),
    .std.resid = residuals(model, type = "standardized")
  )
}

#' @rdname influence.spmodel
#' @method influence spautor
#' @order 2
#' @export
influence.spautor <- influence.splm
