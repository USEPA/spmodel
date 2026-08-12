#' @rdname print.spmodel
#' @method print splmRF
#' @order 13
#' @export
#'
#' @examples
#' \donttest{
#' sulfate$var <- rnorm(NROW(sulfate)) # add noise variable
#' sprfmod <- splmRF(sulfate ~ var, data = sulfate, spcov_type = "exponential")
#' print(sprfmod)
#' }
print.splmRF <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("ranger:\n")
  print(x$ranger)

  cat("\nsplm on ranger residuals:\n")
  print(x$splm, digits = digits, ...)

  invisible(x)
}

#' @rdname print.spmodel
#' @method print summary.splmRF
#' @order 14
#' @export
print.summary.splmRF <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"), ...) {
  cat("ranger:\n")
  print(x$ranger)

  cat("\nsplm on ranger residuals:\n")
  print(x$splm, digits = digits, signif.stars = signif.stars, ...)

  invisible(x)
}

#' @rdname print.spmodel
#' @method print spautorRF
#' @order 15
#' @export
#'
#' @examples
#' \donttest{
#' sprfmod <- spautorRF(log_trend ~ stock, data = seal, spcov_type = "car")
#' print(sprfmod)
#' }
print.spautorRF <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("ranger:\n")
  print(x$ranger)

  cat("\nspautor on ranger residuals:\n")
  print(x$spautor, digits = digits, ...)

  invisible(x)
}

#' @rdname print.spmodel
#' @method print summary.spautorRF
#' @order 16
#' @export
print.summary.spautorRF <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     signif.stars = getOption("show.signif.stars"), ...) {
  cat("ranger:\n")
  print(x$ranger)

  cat("\nspautor on ranger residuals:\n")
  print(x$spautor, digits = digits, signif.stars = signif.stars, ...)

  invisible(x)
}
