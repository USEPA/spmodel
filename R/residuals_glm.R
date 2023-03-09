#' @rdname residuals.spmodel
#' @method residuals spglm
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
#' @export
resid.spglm <- residuals.spglm

#' @rdname residuals.spmodel
#' @method residuals spgautor
#' @export
residuals.spgautor <- residuals.spglm

#' @rdname residuals.spmodel
#' @method resid spgautor
#' @export
resid.spgautor <- residuals.spgautor

#' @rdname rstandard.spmodel
#' @method rstandard spglm
#' @export
rstandard.spglm <- function(model, ...) {
  residuals.spglm(model, type = "standardized")
}

#' @rdname rstandard.spmodel
#' @method rstandard spgautor
#' @export
rstandard.spgautor <- rstandard.spglm
