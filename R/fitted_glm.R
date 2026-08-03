#' @rdname fitted.spmodel
#' @method fitted spglm
#' @export
fitted.spglm <- function(object, type = "response", ...) {
  # these fitted values were already computed and stored on the model object
  # during fitting, so this is a lookup rather than a fresh computation
  # unlike the linear (splm) case, glm models also expose "link"-scale fitted
  # values (i.e., before the inverse link is applied to reach the response scale)
  if (type == "link") {
    fitted_val <- object$fitted$link
  } else if (type == "response") {
    fitted_val <- object$fitted$response
  } else if (type == "spcov") {
    fitted_val <- object$fitted$spcov
  } else if (type == "randcov") {
    fitted_val <- object$fitted$randcov
  } else {
    stop("Invalid type argument. The type argument must be \"response\",  \"link\", \"spcov\", or \"randcov\".", call. = FALSE)
  }
  fitted_val
}

#' @rdname fitted.spmodel
#' @method fitted.values spglm
#' @export
fitted.values.spglm <- fitted.spglm

#' @rdname fitted.spmodel
#' @method fitted spgautor
#' @export
fitted.spgautor <- fitted.spglm

#' @rdname fitted.spmodel
#' @method fitted.values spgautor
#' @export
fitted.values.spgautor <- fitted.spgautor
