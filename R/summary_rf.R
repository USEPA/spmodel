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
  # see print.splmRF() for why both element names are checked
  splm_out <- if (!is.null(object$splm)) object$splm else object$spautor
  summary_list <- list(ranger = object$ranger, splm = summary(splm_out))
  structure(summary_list, class = "summary.splmRF")
}
