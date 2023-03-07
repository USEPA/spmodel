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

#' @export
resid.spglm <- residuals.spglm

#' @export
residuals.spgautor <- residuals.spglm

#' @export
resid.spgautor <- residuals.spgautor

#' @export
rstandard.spglm <- function(model, ...) {
  residuals.spglm(model, type = "standardized")
}

#' @export
rstandard.spgautor <- rstandard.spglm
