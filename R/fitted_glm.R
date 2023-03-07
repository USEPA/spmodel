#' @export
fitted.spglm <- function(object, type = "link", ...) {
  if (type == "link") {
    fitted_val <- object$fitted$link
  } else if (type == "response") {
    fitted_val <- object$fitted$response
  } else {
    stop("Invalid type argument. The type argument must be \"link\" or \"response\".", call. = FALSE)
  }
  fitted_val
}

#' @export
fitted.spgautor <- fitted.spglm
