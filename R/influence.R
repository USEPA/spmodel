#' Regression diagnostics
#'
#' @description Provides basic quantities which are used in forming
#'   a wide variety of diagnostics for checking the quality of fitted model objects.
#'
#' @param model A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details This function calls [residuals.spmod()], [hatvalues.spmod()],
#'   and [cooks.distance.spmod()] and puts the results into a tibble. It is
#'   primarily used when calling [augment.spmod()].
#'
#' @return A tibble with residuals (\code{.resid}), leverage values (\code{.hat}),
#'   cook's distance (\code{.cooksd}), and standardized residuals (\code{.std.resid}).
#'
#' @method influence spmod
#' @export
#'
#' @seealso [augment.spmod()] [cooks.distance.spmod()] [hatvalues.spmod()] [residuals.spmod()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' influence(spmod)
influence.spmod <- function(model, ...) {
  tibble::tibble( # used to be data.frame
    .resid = residuals(model),
    .hat = hatvalues(model),
    .cooksd = cooks.distance(model),
    .std.resid = residuals(model, type = "standardized") # ,
    # .sigma = abs(model$model$y - loocv(model, cv_fitted = TRUE)$cv_fitted)
  )
}
