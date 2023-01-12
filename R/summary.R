#' Summarize a fitted model object
#'
#' @description Summarize a fitted model object.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details \code{summary()} creates a summary of a fitted model object
#'   intended to be printed using \code{print()}. This summary contains
#'   useful information like the original function call, residuals,
#'   a coefficients table, a pseudo r-squared, and estimated covariance
#'   parameters.
#'
#' @return A list with several fitted model quantities used to create
#'   informative summaries when printing.
#'
#' @method summary splm
#' @export
#'
#' @seealso [print.summary.splm()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' summary(spmod)
summary.splm <- function(object, ...) {
  summary_coefficients_fixed <- data.frame(
    estimates = coef(object, type = "fixed"),
    Std_Error = sqrt(diag(vcov(object, type = "fixed")))
  )

  summary_coefficients_fixed$z_value <- summary_coefficients_fixed$estimates / summary_coefficients_fixed$Std_Error
  summary_coefficients_fixed$p <- 2 * (1 - pnorm(abs(summary_coefficients_fixed$z_value)))

  spcov_params_val <- coef(object, type = "spcov")
  randcov_params_val <- coef(object, type = "randcov")
  coefficients <- list(fixed = summary_coefficients_fixed, spcov = spcov_params_val, randcov = randcov_params_val)
  summary_list <- list(
    call = object$call,
    terms = object$terms,
    residuals = object$residuals,
    coefficients = coefficients,
    pseudoR2 = object$pseudoR2,
    vcov = object$vcov,
    is_known = object$is_known,
    anisotropy = object$anisotropy
  )
  new_summary_list <- structure(summary_list, class = paste("summary", class(object), sep = "."))
  new_summary_list
}

#' @rdname summary.splm
#' @method summary spautor
#' @export
summary.spautor <- summary.splm
