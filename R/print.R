#' Print values
#'
#' @description Print fitted model objects and summaries.
#'
#' @param x A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()] or output from
#'   \code{summary(x)} or or \code{anova(x)}.
#' @param digits The number of significant digits to use when printing.
#' @param signif.stars Logical. If \code{TRUE}, significance stars are printed for each coefficient
#' @param ... Other arguments passed to or from other methods.
#'
#' @return Printed fitted model objects and summaries with formatting.
#'
#' @name print.spmodel
#' @method print splm
#' @order 1
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

  spcoef <- select_print_spcoef_pointref(coef(x, type = "spcov"), x$anisotropy)

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
#' @order 2
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

  spcoef <- select_print_spcoef_areal(coef(x, type = "spcov"), x$is_known$spcov)

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
#' @order 3
#' @export
print.summary.splm <- function(x,
                               digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  print_residual_summary(x$residuals, digits, field = "response")

  # pasting the fixed coefficient summary
  cat("\nCoefficients (fixed):\n")
  coefs_fixed <- x$coefficients$fixed
  if ("df" %in% colnames(coefs_fixed)) {
    colnames(coefs_fixed) <- c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")
  } else {
    colnames(coefs_fixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }
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

  if (length(x$coefficients$randcov)) {
    cat("\nCoefficients (random effects):\n")
    print(x$coefficients$randcov, digits = digits)
    cat("\n")
  }

  invisible(x)
}

#' @rdname print.spmodel
#' @method print summary.spautor
#' @order 4
#' @export
print.summary.spautor <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  print_residual_summary(x$residuals, digits, field = "response")

  # pasting the fixed coefficient summary
  cat("\nCoefficients (fixed):\n")
  coefs_fixed <- x$coefficients$fixed
  if ("df" %in% colnames(coefs_fixed)) {
    colnames(coefs_fixed) <- c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")
  } else {
    colnames(coefs_fixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }
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

  if (length(x$coefficients$randcov)) {
    cat("\nCoefficients (random effects):\n")
    print(x$coefficients$randcov, digits = digits)
    cat("\n")
  }

  invisible(x)
}

#' @rdname print.spmodel
#' @method print anova.splm
#' @order 5
#' @export
print.anova.splm <- function(x, digits = max(getOption("digits") - 2L, 3L),
                             signif.stars = getOption("show.signif.stars"), ...) {
  cat(attr(x, "heading")[1])
  cat("\n")
  cat(attr(x, "heading")[2])
  cat("\n")
  # a Pr(>Chi2)/Pr(>F) column is only present when test = TRUE (Pr(>F)
  # specifically when ddf = "satterthwaite" was used -- see anova.splm())
  if ("Pr(>Chi2)" %in% colnames(x) || "Pr(>F)" %in% colnames(x)) {
    P.values <- TRUE
    has.Pvalue <- TRUE
  } else {
    P.values <- FALSE
    has.Pvalue <- FALSE
  }
  if ("NumDF" %in% colnames(x)) {
    # NumDF (always a whole number) is otherwise grouped with DenDF under
    # printCoefmat()'s default cs.ind, which rounds both to a shared decimal
    # precision and puts spurious trailing zeros (e.g. "1.000") on NumDF;
    # pointing cs.ind/tst.ind at DenDF/F value alone leaves NumDF to fall
    # through to plain format(), printing as an integer
    cs.ind <- which(colnames(x) == "DenDF")
    tst.ind <- which(colnames(x) == "F value")
    printCoefmat(x, digits = digits, signif.stars = signif.stars, P.values = P.values, has.Pvalue = has.Pvalue, cs.ind = cs.ind, tst.ind = tst.ind, ...)
  } else {
    printCoefmat(x, digits = digits, signif.stars = signif.stars, P.values = P.values, has.Pvalue = has.Pvalue, ...)
  }
}

#' @rdname print.spmodel
#' @method print anova.spautor
#' @order 6
#' @export
print.anova.spautor <- print.anova.splm
