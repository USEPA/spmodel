#' Title
#'
#' @param object filler
#' @param newdata filler
#' @param ... filler
#'
#' @rdname predict.spmodel
#' @method predict decorrelate
#'
#' @return filler
#' @export
predict.decorrelate <- function(object, newdata, local, ...) {

  if (missing(local)) local <- object$decorrelate_data$local
  tnewdata <- decorrelate_newdata(object$decorrelate_data, newdata, local, ...)
  tpreds <- predict_decorrelate_algorithm(object$fit, tnewdata, object$algorithm)
  preds <- recorrelate_newdata(tpreds, tnewdata)
  preds
}
