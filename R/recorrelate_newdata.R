#' Recorrelate Machine Learning Predictions
#'
#' @description Recorrelate machine learning predictions according to a
#'   spatial decorrelation transformation.
#'
#'
#' @param object A [decorrelate_newdata()] object.
#' @param ty_newdata Predictions from the machine learning model trained on the
#'   spatially decorrelated data and applied to spatially decorrelated newdata.
#'
#' @return A vector of predictions
#' @export
#'
#' @examples
#' params <- spcov_params("exponential", de = 1, ie = 0.2, range = 1e5)
#' decorr <- decorrelate_data(log_cond ~ temp, data = lake, spcov_params = params)
#' fit <- ranger::ranger(x = decorr$tX, y = decorr$ty)
#' decorr_newdata <- decorrelate_newdata(decorr, newdata = lake_preds)
#' rfpreds <- predict(fit, data = decorr_newdata$tX_newdata)$predictions
#' recorrelate_newdata(decorr_newdata, rfpreds)
recorrelate_newdata <- function(object, ty_newdata) {

  if (!inherits(object, "decorrelate_newdata")) {
    stop("object must have class \"decorrelate_newdata\".", call. = FALSE)
  }

  # inverts the response-scale part of the spatial decorrelation transform:
  # object$yscale/yoffset are the per-observation conditional standard
  # deviation/mean computed by get_decorrelate_newdata() when object was built
  output <- object$yscale * ty_newdata + object$yoffset
  # object$y was built on the offset-subtracted scale (get_data_object_splm()
  # subtracts any formula offset() term before decorrelation), so newdata's
  # own offset must be added back here to return to the response scale --
  # the same "subtract at the start, add back at the end" pattern
  # conditional.splm()/predict.splm() use
  # confusingly, it is important that yoffset is the part added back
  # in the recorrelation while offset is the standard offset term
  if (!is.null(object$offset)) {
    output <- output + object$offset
  }
  output

}
