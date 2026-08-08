#' @rdname print.spmodel
#' @method print spglm
#' @order 7
#' @export
print.spglm <- function(x, digits = max(3L, getOption("digits") - 3L),
                        ...) {
  cat("\nCall:\n", paste(deparse(x$call),
    sep = "\n",
    collapse = "\n"
  ), "\n\n", sep = "")

  cat("\n")

  cat("Coefficients (fixed):\n")
  print.default(format(coef(x, type = "fixed"), digits = digits),
    print.gap = 2L,
    quote = FALSE
  )

  cat("\n")

  spcoef <- select_print_spcoef_pointref(coef(x, type = "spcov"), x$anisotropy)

  cat(paste("\nCoefficients (", class(coef(x, type = "spcov")), " spatial covariance):\n", sep = ""))
  print.default(format(spcoef, digits = digits),
    print.gap = 2L,
    quote = FALSE
  )

  cat("\n")

  # dispersion is the extra family-specific scale/shape parameter (e.g.
  # overdispersion) that generalized linear models need beyond the mean structure
  print_dispersion_block(coef(x, type = "dispersion"), digits, style = "raw")

  if (length(coef(x, type = "randcov"))) {
    cat("Coefficients (random effects):\n")
    print.default(format(coef(x, type = "randcov"), digits = digits),
      print.gap = 2L,
      quote = FALSE
    )

    cat("\n")
  }
  invisible(x)
}

#' @rdname print.spmodel
#' @method print spgautor
#' @order 8
#' @export
print.spgautor <- function(x, digits = max(3L, getOption("digits") - 3L),
                           ...) {
  cat("\nCall:\n", paste(deparse(x$call),
    sep = "\n",
    collapse = "\n"
  ), "\n\n", sep = "")

  cat("\n")

  cat("Coefficients (fixed):\n")
  print.default(format(coef(x, type = "fixed"), digits = digits),
    print.gap = 2L,
    quote = FALSE
  )

  cat("\n")

  spcoef <- select_print_spcoef_areal(coef(x, type = "spcov"), x$is_known$spcov)

  cat(paste("\nCoefficients (", class(coef(x, type = "spcov")), " spatial covariance):\n", sep = ""))
  print.default(format(spcoef, digits = digits),
    print.gap = 2L,
    quote = FALSE
  )

  cat("\n")

  print_dispersion_block(coef(x, type = "dispersion"), digits, style = "raw")

  if (length(coef(x, type = "randcov"))) {
    cat("Coefficients (random effects):\n")
    print.default(format(coef(x, type = "randcov"), digits = digits),
      print.gap = 2L,
      quote = FALSE
    )

    cat("\n")
  }
  invisible(x)
}

#' @rdname print.spmodel
#' @method print summary.spglm
#' @order 9
#' @export
print.summary.spglm <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  print_residual_summary(x$residuals, digits, field = "deviance")

  # pasting the fixed coefficient summary
  cat("\nCoefficients (fixed):\n")
  coefs_fixed <- x$coefficients$fixed
  colnames(coefs_fixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  printCoefmat(coefs_fixed, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  # pasting the generalized r squared
  if (x$pseudoR2 != 0) {
    cat("\nPseudo R-squared: ")
    cat(formatC(x$pseudoR2, digits = digits))
    cat("\n")
  }

  # pasting the covariance coefficient summary
  spcoef <- select_print_spcoef_pointref(x$coefficients$spcov, x$anisotropy)

  cat(paste("\nCoefficients (", class(x$coefficients$spcov), " spatial covariance):\n", sep = ""))
  print(spcoef, digits = digits)

  print_dispersion_block(x$coefficients$dispersion, digits, style = "summary")

  if (length(x$coefficients$randcov)) {
    cat("\nCoefficients (random effects):\n")
    print(x$coefficients$randcov, digits = digits)
    cat("\n")
  }

  invisible(x)
}

#' @rdname print.spmodel
#' @method print summary.spgautor
#' @order 10
#' @export
print.summary.spgautor <- function(x,
                                   digits = max(3L, getOption("digits") - 3L),
                                   signif.stars = getOption("show.signif.stars"),
                                   ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  print_residual_summary(x$residuals, digits, field = "deviance")

  # pasting the fixed coefficient summary
  cat("\nCoefficients (fixed):\n")
  coefs_fixed <- x$coefficients$fixed
  colnames(coefs_fixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  printCoefmat(coefs_fixed, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  # pasting the generalized r squared
  if (x$pseudoR2 != 0) {
    cat("\nPseudo R-squared: ")
    cat(formatC(x$pseudoR2, digits = digits))
    cat("\n")
  }

  # pasting the covariance coefficient summary
  spcoef <- select_print_spcoef_areal(x$coefficients$spcov, x$is_known$spcov)

  cat(paste("\nCoefficients (", class(x$coefficients$spcov), " spatial covariance):\n", sep = ""))
  print(spcoef, digits = digits)

  print_dispersion_block(x$coefficients$dispersion, digits, style = "summary")

  if (length(x$coefficients$randcov)) {
    cat("\nCoefficients (random effects):\n")
    print(x$coefficients$randcov, digits = digits)
    cat("\n")
  }

  invisible(x)
}

#' @rdname print.spmodel
#' @method print anova.spglm
#' @order 11
#' @export
print.anova.spglm <- print.anova.splm

#' @rdname print.spmodel
#' @method print anova.spgautor
#' @order 12
#' @export
print.anova.spgautor <- print.anova.spautor
