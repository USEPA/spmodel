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

  output <- object$yscale * ty_newdata + object$yoffset
  output

}
