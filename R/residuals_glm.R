#' @rdname residuals.spmodel
#' @method residuals spglm
#' @order 7
#' @export
residuals.spglm <- function(object, type = "deviance", ...) {
  if (type == "deviance") {
    return(object$residuals$deviance)
  } else if (type == "response") {
    return(object$residuals$response)
  } else if (type == "pearson") {
    return(object$residuals$pearson)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("residuals must be deviance or response or pearson or standardized")
  }
}

#' @rdname residuals.spmodel
#' @method resid spglm
#' @order 8
#' @export
resid.spglm <- residuals.spglm

#' @rdname residuals.spmodel
#' @method residuals spgautor
#' @order 10
#' @export
residuals.spgautor <- residuals.spglm

#' @rdname residuals.spmodel
#' @method resid spgautor
#' @order 11
#' @export
resid.spgautor <- residuals.spgautor

#' @rdname residuals.spmodel
#' @method rstandard spglm
#' @order 9
#' @export
rstandard.spglm <- function(model, ...) {
  residuals.spglm(model, type = "standardized")
}

#' @rdname residuals.spmodel
#' @method rstandard spgautor
#' @order 12
#' @export
rstandard.spgautor <- rstandard.spglm
