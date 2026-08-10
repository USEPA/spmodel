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
  # splmRF() always names the residual model element "splm"; spautorRF()
  # can also produce a "splmRF"-classed object (spautor() delegates to splm()
  # internally when spcov_type is "none"/"ie"), but keeps calling that element
  # "spautor" for consistency with its own output naming -- fall back to it
  # so this method works for either constructor
  splm_out <- if (!is.null(x$splm)) x$splm else x$spautor

  cat("ranger:\n")
  print(x$ranger)

  cat("\nsplm on ranger residuals:\n")
  print(splm_out, digits = digits, ...)

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
