# shared building blocks for print.R / print_glm.R
#' Select which spatial covariance coefficients to print for point-referenced models
#'
#' @param spcoef A \code{spcov_params} object (i.e. the output of
#'   \code{coef(x, type = "spcov")} or \code{x$coefficients$spcov})
#' @param anisotropy Whether anisotropy was modeled
#'
#' @return \code{spcoef}, with \code{rotate}/\code{scale} dropped unless
#'   \code{anisotropy} is \code{TRUE}, and collapsed to just \code{ie} when
#'   the covariance type is \code{"none"} or \code{"ie"}
#'
#' @details The \code{"none"}/\code{"ie"} check is captured from \code{spcoef}
#'   *before* any subsetting, because \code{[} subsetting drops a plain
#'   vector's class attribute -- checking \code{inherits()} after the
#'   \code{rotate}/\code{scale} subset would silently always be \code{FALSE}.
#'
#' @noRd
select_print_spcoef_pointref <- function(spcoef, anisotropy) {
  is_none_ie <- inherits(spcoef, c("none", "ie"))
  if (!anisotropy) {
    spcoef <- spcoef[-which(names(spcoef) %in% c("rotate", "scale"))]
  }
  if (is_none_ie) {
    spcoef <- spcoef["ie"]
  }
  spcoef
}

#' Select which spatial covariance coefficients to print for areal models
#'
#' @param spcoef A \code{spcov_params} object (i.e. the output of
#'   \code{coef(x, type = "spcov")} or \code{x$coefficients$spcov})
#' @param is_known_spcov The \code{is_known$spcov} element of the fitted
#'   model object (or its summary), a named logical vector
#'
#' @return \code{spcoef}, subset to \code{de}/\code{range} plus whichever of
#'   \code{ie}/\code{extra} were not fixed at a known value of zero
#'
#' @noRd
select_print_spcoef_areal <- function(spcoef, is_known_spcov) {
  no_ie <- spcoef[["ie"]] == 0 && is_known_spcov[["ie"]]
  no_extra <- spcoef[["extra"]] == 0 && is_known_spcov[["extra"]]
  if (no_ie && no_extra) {
    spcoef[c("de", "range")]
  } else if (no_ie && !no_extra) {
    spcoef[c("de", "range", "extra")]
  } else if (!no_ie && no_extra) {
    spcoef[c("de", "ie", "range")]
  } else {
    spcoef[c("de", "ie", "range", "extra")]
  }
}

#' Print the dispersion coefficient block for GLM-type model objects
#'
#' @param dispersion_val A \code{family}-classed dispersion coefficient
#'   object (i.e. \code{coef(x, type = "dispersion")} or
#'   \code{x$coefficients$dispersion})
#' @param digits The number of significant digits to use
#' @param style Either \code{"raw"} (fitted model object print methods,
#'   which use \code{print.default(..., print.gap = 2L, quote = FALSE)} and a
#'   trailing blank line) or \code{"summary"} (summary print methods, which
#'   use plain \code{print()} and no trailing blank line) -- these are not
#'   unified into one rendering because they do not currently produce
#'   byte-identical output (confirmed: \code{print.gap = 2L} adds an extra
#'   trailing space per column versus the default gap), matching the same
#'   raw-vs-summary asymmetry base R's own \code{print.lm()}/
#'   \code{print.summary.lm()} have.
#'
#' @noRd
print_dispersion_block <- function(dispersion_val, digits, style = c("raw", "summary")) {
  style <- match.arg(style)
  cat(paste("\nCoefficients (Dispersion for ", class(dispersion_val), " family):\n", sep = ""))
  if (style == "raw") {
    print.default(format(unclass(dispersion_val), digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
    cat("\n")
  } else {
    print(unclass(dispersion_val), digits = digits)
  }
}

#' Print the residual quantile summary block for summary print methods
#'
#' @param residuals_list The \code{residuals} element of a fitted model
#'   object or its summary (a list with elements including \code{response}
#'   and, for GLM-type objects, \code{deviance})
#' @param digits The number of significant digits to use
#' @param field Either \code{"response"} (Gaussian models, header
#'   \code{"Residuals:"}) or \code{"deviance"} (GLM-type models, header
#'   \code{"Deviance Residuals:"})
#'
#' @noRd
print_residual_summary <- function(residuals_list, digits, field = c("response", "deviance")) {
  field <- match.arg(field)
  header <- if (field == "deviance") "Deviance Residuals:" else "Residuals:"
  resid_vals <- residuals_list[[field]]
  cat(paste0("\n", header, "\n"))
  resQ <- c(
    min(resid_vals), quantile(resid_vals, p = c(0.25, 0.5, 0.75), na.rm = TRUE),
    max(resid_vals)
  )
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)
}
