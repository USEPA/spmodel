#' @rdname influence.spmodel
#' @method influence spglm
#' @order 3
#' @export
influence.spglm <- function(model, ...) {
  tibble::tibble(
    # allow ... so type.residuals can be passed to residuals()
    .resid = residuals(model, ...),
    .hat = hatvalues(model),
    .cooksd = cooks.distance(model),
    # standardized residuals always use their own default type regardless of
    # any type passed above, since standardization needs a specific residual definition
    .std.resid = residuals(model, type = "standardized")
  )
}

#' @rdname influence.spmodel
#' @method influence spgautor
#' @order 4
#' @export
influence.spgautor <- influence.spglm
