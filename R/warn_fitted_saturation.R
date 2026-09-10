#' Warn when nearly all fitted binomial probabilities are saturated near 0 or 1
#'
#' Spatial structure can make binomial \code{spglm()}/\code{spgautor()} fits
#' separate far more readily than an ordinary (non-spatial) logistic
#' regression, even with well-behaved covariates. Deliberately not
#' \code{glm.fit}'s own \code{any(mu > 1 - eps) | any(mu < eps)} convention at a
#' loose tolerance: for heavily-imbalanced data, that convention already fires
#' on a perfectly reasonable non-spatial fit, so it does not discriminate a
#' pathological fit from a fine one. Requiring nearly all fitted probabilities
#' to be saturated at a tight tolerance tracks the standard-error blowup that
#' accompanies this kind of separation far more specifically.
#'
#' @noRd
warn_fitted_saturation <- function(fitted_response, family) {
  if (family != "binomial") {
    return(invisible())
  }
  tol <- 1e-6
  saturated <- fitted_response < tol | fitted_response > 1 - tol
  if (mean(saturated) >= 0.99) {
    warning(
      "Nearly all fitted probabilities are numerically 0 or 1. Perfect separation detected.",
      call. = FALSE
    )
  }
  invisible()
}
