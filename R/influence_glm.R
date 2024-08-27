#' @rdname influence.spmodel
#' @method influence spglm
#' @order 3
#' @export
influence.spglm <- function(model, ...) {
  tibble::tibble( # used to be data.frame
    # allow ... so type.residuals can be passed to residuals()
    .resid = residuals(model, ...),
    .hat = hatvalues(model),
    .cooksd = cooks.distance(model),
    .std.resid = residuals(model, type = "standardized") # ,
    # .sigma = abs(model$model$y - loocv(model, cv_fitted = TRUE)$cv_fitted)
  )
}

#' @rdname influence.spmodel
#' @method influence spgautor
#' @order 4
#' @export
influence.spgautor <- influence.spglm
