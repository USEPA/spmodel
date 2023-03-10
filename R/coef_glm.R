#' @rdname coef.spmodel
#' @method coef spglm
#' @export
coef.spglm <- function(object, type = "fixed", ...) {
  if (type == "fixed") {
    return(object$coefficients$fixed)
  } else if (type == "spcov") {
    spcov_coef <- object$coefficients$spcov
    return(spcov_coef)
  } else if (type == "dispersion") {
    return(object$coefficients$dispersion)
  } else if (type == "randcov") {
    return(object$coefficients$randcov)
  } else {
    stop("Invalid type argument. The type argument must be \"fixed\", \"spcov\", \"dispersion\", or \"randcov\".", call. = FALSE)
  }
}

#' @rdname coef.spmodel
#' @method coefficients spglm
#' @export
coefficients.spglm <- coef.spglm

#' @rdname coef.spmodel
#' @method coef spgautor
#' @export
coef.spgautor <- coef.spglm

#' @rdname coef.spmodel
#' @method coefficients spgautor
#' @export
coefficients.spgautor <- coef.spgautor
