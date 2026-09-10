#' Summarize a fitted model object
#'
#' @description Summarize a fitted model object.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], [spgautor()], [splmRF()], or [spautorRF()].
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
#' @name summary.spmodel
#' @method summary splm
#' @order 1
#' @export
#'
#' @seealso [print.spmodel()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' summary(spmod)
summary.splm <- function(object, ...) {
  # standard errors come from the square root of the diagonal of the fixed
  # effect covariance matrix (off-diagonal covariances are ignored here,
  # as is standard when reporting per-coefficient standard errors)
  summary_coefficients_fixed <- data.frame(
    estimates = coef(object, type = "fixed"),
    Std_Error = sqrt(diag(vcov(object, type = "fixed")))
  )

  if (!is.null(object$ddf)) {
    # denominator degrees of freedom available (see splm()/spautor()'s ddf
    # argument) -- t-test each fixed effect using them instead of the
    # asymptotic z-test below, matching lmerTest's Satterthwaite summary()
    summary_coefficients_fixed$df <- object$ddf[rownames(summary_coefficients_fixed)]
    summary_coefficients_fixed$t_value <- summary_coefficients_fixed$estimates / summary_coefficients_fixed$Std_Error
    summary_coefficients_fixed$p <- 2 * pt(abs(summary_coefficients_fixed$t_value), summary_coefficients_fixed$df, lower.tail = FALSE)
  } else {
    # Wald z-test for each fixed effect: estimate / se ~ N(0, 1) under the null
    # that the true coefficient is zero, so the two-sided p-value uses the
    # standard normal (infinite t df) distribution
    summary_coefficients_fixed$z_value <- summary_coefficients_fixed$estimates / summary_coefficients_fixed$Std_Error
    summary_coefficients_fixed$p <- 2 * (1 - pnorm(abs(summary_coefficients_fixed$z_value)))
  }

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
  # tag the class as "summary.<original class>" (e.g. summary.splm) so
  # print.spmodel() can dispatch on the fitted model's type
  new_summary_list <- structure(summary_list, class = paste("summary", class(object), sep = "."))
  new_summary_list
}

#' @rdname summary.spmodel
#' @method summary spautor
#' @order 2
#' @export
summary.spautor <- summary.splm
