#' @rdname predict.spmodel
#' @method predict decorrelate
#' @order 13
#' @export
#' @examples
#' decorr <- decorrelate(log_cond ~ temp, data = lake, spcov_type = "exponential")
#' predict(decorr, newdata = lake_preds)
predict.decorrelate <- function(object, newdata, local, ...) {


  if (missing(newdata) || is.null(newdata)) newdata <- object$newdata # always exists as object or NULL
  # error if newdata missing from arguments and object
  if (is.null(newdata) || NROW(newdata) == 0) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }
  if (missing(local) || is.null(local)) local <- object$decorrelate_data$local
  tnewdata <- decorrelate_newdata(object$decorrelate_data, newdata, local, ...)
  tpreds <- predict_decorrelate_algorithm(object$fit, tnewdata, object$algorithm)
  preds <- recorrelate_newdata(tnewdata, tpreds)
  preds
}

#' @rdname predict.spmodel
#' @method predict decorrelate_list
#' @order 14
#' @export
predict.decorrelate_list <- function(object, newdata, local, ...) {

  if (missing(newdata)) newdata <- NULL
  if (missing(local)) local <- NULL
  preds <- lapply(object, function(x) predict(x, newdata, local, ...))
  names(preds) <- names(object)
  preds
}
