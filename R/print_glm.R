#' @method print spglm
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

  spcoef <- coef(x, type = "spcov")

  if (!x$anisotropy) {
    spcoef <- spcoef[-which(names(spcoef) %in% c("rotate", "scale"))]
  }
  if (inherits(coef(x, type = "spcov"), "none")) {
    spcoef <- spcoef["ie"]
  }

  cat(paste("\nCoefficients (", class(coef(x, type = "spcov")), " spatial covariance):\n", sep = ""))
  print.default(format(spcoef, digits = digits),
                print.gap = 2L,
                quote = FALSE
  )

  cat("\n")

  cat(paste("\nCoefficients (Dispersion for ", class(coef(x, type = "dispersion")), " family):\n", sep = ""))
  print.default(format(unclass(coef(x, type = "dispersion")), digits = digits),
                print.gap = 2L,
                quote = FALSE
  )

  cat("\n")

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

  spcoef <- coef(x, type = "spcov")

  no_ie <- spcoef[["ie"]] == 0 && x$is_known$spcov[["ie"]]
  no_extra <- spcoef[["extra"]] == 0 && x$is_known$spcov[["extra"]]
  if (no_ie && no_extra) {
    spcoef <- spcoef[c("de", "range")]
  } else if (no_ie && !no_extra) {
    spcoef <- spcoef[c("de", "range", "extra")]
  } else if (!no_ie && no_extra) {
    spcoef <- spcoef[c("de", "ie", "range")]
  } else {
    spcoef <- spcoef[c("de", "ie", "range", "extra")]
  }

  cat(paste("\nCoefficients (", class(coef(x, type = "spcov")), " spatial covariance):\n", sep = ""))
  print.default(format(spcoef, digits = digits),
                print.gap = 2L,
                quote = FALSE
  )

  cat("\n")

  cat(paste("\nCoefficients (Dispersion for ", class(coef(x, type = "dispersion")), " family):\n", sep = ""))
  print.default(format(unclass(coef(x, type = "dispersion")), digits = digits),
                print.gap = 2L,
                quote = FALSE
  )

  cat("\n")

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

#' @export
print.summary.spglm <- function(x,
                               digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  cat("\nResiduals:\n")
  resQ <- c(
    min(x$residuals$response), quantile(x$residuals$response, p = c(0.25, 0.5, 0.75), na.rm = TRUE),
    max(x$residuals$response)
  )
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)

  # pasting the fixed coefficient summary
  cat("\nCoefficients (fixed):\n")
  coefs_fixed <- x$coefficients$fixed
  # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  colnames(coefs_fixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  printCoefmat(coefs_fixed, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  # pasting the generalized r squared
  if (x$pseudoR2 != 0) {
    cat("\nPseudo R-squared: ")
    cat(formatC(x$pseudoR2, digits = digits))
    cat("\n")
  }

  # pasting the covariance coefficient summary
  spcoef <- x$coefficients$spcov


  if (!x$anisotropy) {
    spcoef <- spcoef[-which(names(spcoef) %in% c("rotate", "scale"))]
  }
  if (inherits(x$coefficients$spcov, "none")) {
    spcoef <- spcoef["ie"]
  }


  cat(paste("\nCoefficients (", class(x$coefficients$spcov), " spatial covariance):\n", sep = ""))
  print(spcoef, digits = digits)

  cat(paste("\nCoefficients (Dispersion for ", class(x$coefficients$dispersion), " family):\n", sep = ""))
  print(unclass(x$coefficients$dispersion), digits = digits)

  if (length(x$coefficients$randcov)) {
    cat("\nCoefficients (random effects):\n")
    print(x$coefficients$randcov, digits = digits)
    cat("\n")
  }

  invisible(x)
}

#' @export
print.summary.spgautor <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  cat("\nResiduals:\n")
  resQ <- c(
    min(x$residuals$response), quantile(x$residuals$response, p = c(0.25, 0.5, 0.75), na.rm = TRUE),
    max(x$residuals$response)
  )
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)

  # pasting the fixed coefficient summary
  cat("\nCoefficients (fixed):\n")
  coefs_fixed <- x$coefficients$fixed
  # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  colnames(coefs_fixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  printCoefmat(coefs_fixed, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  # pasting the generalized r squared
  if (x$pseudoR2 != 0) {
    cat("\nPseudo R-squared: ")
    cat(formatC(x$pseudoR2, digits = digits))
    cat("\n")
  }

  # pasting the covariance coefficient summary
  spcoef <- x$coefficients$spcov


  no_ie <- spcoef[["ie"]] == 0 && x$is_known$spcov[["ie"]]
  no_extra <- spcoef[["extra"]] == 0 && x$is_known$spcov[["extra"]]
  if (no_ie && no_extra) {
    spcoef <- spcoef[c("de", "range")]
  } else if (no_ie && !no_extra) {
    spcoef <- spcoef[c("de", "range", "extra")]
  } else if (!no_ie && no_extra) {
    spcoef <- spcoef[c("de", "ie", "range")]
  } else {
    spcoef <- spcoef[c("de", "ie", "range", "extra")]
  }


  cat(paste("\nCoefficients (", class(x$coefficients$spcov), " spatial covariance):\n", sep = ""))
  print(spcoef, digits = digits)

  cat(paste("\nCoefficients (Dispersion for ", class(x$coefficients$dispersion), " family):\n", sep = ""))
  print(unclass(x$coefficients$dispersion), digits = digits)

  if (length(x$coefficients$randcov)) {
    cat("\nCoefficients (random effects):\n")
    print(x$coefficients$randcov, digits = digits)
    cat("\n")
  }

  invisible(x)
}

#' @export
print.anova.spglm <- print.anova.splm

#' @export
print.anova.spgautor <- print.anova.spautor
