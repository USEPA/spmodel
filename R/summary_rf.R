#' @rdname summary.spmodel
#' @method summary splmRF
#' @order 5
#' @export
#'
#' @examples
#' \donttest{
#' sprfmod <- splmRF(log_cond ~ temp + precip, data = lake, spcov_type = "exponential")
#' summary(sprfmod)
#' }
summary.splmRF <- function(object, ...) {
  summary_list <- list(ranger = object$ranger, splm = summary(object$splm))
  structure(summary_list, class = "summary.splmRF")
}

#' @rdname summary.spmodel
#' @method summary spautorRF
#' @order 6
#' @export
#'
#' @examples
#' \donttest{
#' sprfmod <- spautorRF(log_trend ~ stock, data = seal, spcov_type = "car")
#' summary(sprfmod)
#' }
summary.spautorRF <- function(object, ...) {
  summary_list <- list(ranger = object$ranger, spautor = summary(object$spautor))
  structure(summary_list, class = "summary.spautorRF")
}
