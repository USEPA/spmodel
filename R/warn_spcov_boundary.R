#' Warn when an ml fit's spatial covariance has collapsed toward the numerical boundary
#'
#' estmethod = "ml" can report an artificially improved likelihood as de + ie -> 0,
#' making likelihood-based comparisons (AIC, BIC, etc.) unreliable near this
#' boundary -- including for spcov_type = "none" and for comparisons across
#' different fixed effects (see the ml-boundary write-up for the full mechanism).
#'
#' (need to update to account for random effects)
#'
#' @noRd
warn_spcov_boundary <- function(spcov_params_val, diagtol) {
  # car/sar never floor ie and have a different (range-based) singularity story;
  # diagtol <= 0 means there's no floor concept to be near
  if (inherits(spcov_params_val, c("car", "sar")) || diagtol <= 0) {
    return(invisible())
  }
  de <- if ("de" %in% names(spcov_params_val)) spcov_params_val[["de"]] else 0
  ie <- spcov_params_val[["ie"]]
  # de AND ie must be jointly small -- spcov_type = "none" (de identical 0)
  # always trips this, as intended
  if (de + ie <= 10 * diagtol) {
    warning(
      "The fitted spatial variance parameters (de and ie) are near a numerical boundary of zero,
      so the likelihood value may be unreliable for model comparisons (e.g., AIC(), AICc(), BIC()).",
      call. = FALSE
    )
  }
  invisible()
}
