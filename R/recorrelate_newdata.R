recorrelate_newdata <- function(ty_newdata, object) {

  if (!inherits(object, "decorrelate_newdata")) {
    stop("object must have class \"decorrelate_newdata\".", call. = FALSE)
  }

  output <- object$yscale * ty_newdata + object$yoffset
  output

}
