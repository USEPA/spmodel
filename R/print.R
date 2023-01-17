#' Print values
#'
#' @description Print fitted model objects and summaries.
#'
#' @param x A fitted model object from [splm()] or [spautor()] or output from
#'   \code{summary(x)} or or \code{anova(x)}.
#' @param digits The number of significant digits to use when printing.
#' @param signif.stars Logical. If \code{TRUE}, significance stars are printed for each coefficient
#' @param ... Other arguments passed to or from other methods.
#'
#' @return Printed fitted model objects and summaries with formatting.
#'
#' @name print.spmodel
#' @method print splm
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' print(spmod)
#' print(summary(spmod))
#' print(anova(spmod))
print.splm <- function(x, digits = max(3L, getOption("digits") - 3L),
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
#' @method print spautor
#' @export
print.spautor <- function(x, digits = max(3L, getOption("digits") - 3L),
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
#' @method print summary.splm
#' @export
print.summary.splm <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  cat("\nResiduals:\n")
  resQ <- c(
    min(x$residuals$raw), quantile(x$residuals$raw, p = c(0.25, 0.5, 0.75), na.rm = TRUE),
    max(x$residuals$raw)
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

  if (length(x$coefficients$randcov)) {
    cat("\nCoefficients (random effects):\n")
    print(x$coefficients$randcov, digits = digits)
    cat("\n")
  }

  invisible(x)
}

#' @rdname print.spmodel
#' @method print summary.spautor
#' @export
print.summary.spautor <- function(x,
                               digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  cat("\nResiduals:\n")
  resQ <- c(
    min(x$residuals$raw), quantile(x$residuals$raw, p = c(0.25, 0.5, 0.75), na.rm = TRUE),
    max(x$residuals$raw)
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

  if (length(x$coefficients$randcov)) {
    cat("\nCoefficients (random effects):\n")
    print(x$coefficients$randcov, digits = digits)
    cat("\n")
  }

  invisible(x)
}

#' @rdname print.spmodel
#' @method print anova.splm
#' @export
print.anova.splm <- function(x, digits = max(getOption("digits") - 2L, 3L),
                              signif.stars = getOption("show.signif.stars"), ...) {
  cat(attr(x, "heading")[1])
  cat("\n")
  cat(attr(x, "heading")[2])
  cat("\n")
  if ("Pr(>Chi2)" %in% colnames(x)) {
    P.values <- TRUE
    has.Pvalue <- TRUE
  } else {
    P.values <- FALSE
    has.Pvalue <- FALSE
  }
  printCoefmat(x, digits = digits, signif.stars = signif.stars, P.values = P.values, has.Pvalue = has.Pvalue, ...)
}

#' @rdname print.spmodel
#' @method print anova.spautor
#' @export
print.anova.spautor <- print.anova.splm
