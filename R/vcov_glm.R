#' @rdname vcov.spmodel
#' @method vcov spglm
#' @order 3
#' @export
vcov.spglm <- function(object, var_correct = TRUE, ...) {
  type <- "fixed"
  if (type == "fixed") {
    if (var_correct) {
      return(object$vcov$fixed$corrected)
    } else {
      return(object$vcov$fixed$uncorrected)
    }
  }
}

#' @rdname vcov.spmodel
#' @method vcov spgautor
#' @order 4
#' @export
vcov.spgautor <- vcov.spglm
