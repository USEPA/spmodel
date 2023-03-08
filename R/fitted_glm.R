#' @export
fitted.spglm <- function(object, type = "response", ...) {
  if (type == "link") {
    fitted_val <- object$fitted$link
  } else if (type == "response") {
    fitted_val <- object$fitted$response
  } else {
    stop("Invalid type argument. The type argument must be \"response\" or  \"link\".", call. = FALSE)
  }
  fitted_val
}

#' @export
fitted.spgautor <- fitted.spglm
