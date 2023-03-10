#' @rdname summary.spmodel
#' @method summary spglm
#' @export
summary.spglm <- function(object, ...) {
  summary_coefficients_fixed <- data.frame(
    estimates = coef(object, type = "fixed"),
    Std_Error = sqrt(diag(vcov(object, type = "fixed")))
  )

  summary_coefficients_fixed$z_value <- summary_coefficients_fixed$estimates / summary_coefficients_fixed$Std_Error
  summary_coefficients_fixed$p <- 2 * (1 - pnorm(abs(summary_coefficients_fixed$z_value)))

  spcov_params_val <- coef(object, type = "spcov")
  dispersion_params_val <- coef(object, type = "dispersion")
  randcov_params_val <- coef(object, type = "randcov")
  coefficients <- list(fixed = summary_coefficients_fixed, spcov = spcov_params_val,
                       dispersion = dispersion_params_val, randcov = randcov_params_val)
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

#' @rdname summary.spmodel
#' @method summary spgautor
#' @export
summary.spgautor <- summary.spglm
