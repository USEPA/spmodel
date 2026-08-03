#' Correct the \code{anisotropy} flag based on fixed rotate/scale initial values
#'
#' @param anisotropy The user-supplied (or default) anisotropy logical
#' @param spcov_initial A \code{spcov_initial} object
#'
#' @return \code{anisotropy}, flipped to \code{FALSE} if \code{rotate} and
#'   \code{scale} are both fixed at their isotropic values (0 and 1), or
#'   flipped to \code{TRUE} if either is fixed at a non-isotropic value or is
#'   being estimated
#'
#' @noRd
get_anisotropy_corrected <- function(anisotropy, spcov_initial) {
  # anisotropy = TRUE was requested, but if the user has fixed rotate = 0 and
  # scale = 1 (the isotropic values) as known, there is nothing to estimate,
  # so silently fall back to the isotropic (anisotropy = FALSE) model
  if (anisotropy) {
    if (all(c("rotate", "scale") %in% names(spcov_initial$initial))) {
      is_rotate_zero <- (!is.na(spcov_initial$initial[["rotate"]])) && spcov_initial$initial[["rotate"]] == 0
      is_rotate_known <- spcov_initial$is_known[["rotate"]]
      is_scale_one <- (!is.na(spcov_initial$initial[["scale"]])) && spcov_initial$initial[["scale"]] == 1
      is_scale_known <- spcov_initial$is_known[["scale"]]
      if (is_rotate_zero && is_rotate_known && is_scale_one && is_scale_known) {
        anisotropy <- FALSE
      }
    }
  } else {
    # anisotropy = FALSE was requested, but if the user supplied a rotate or
    # scale initial value that is either not fixed or not at its isotropic
    # value, anisotropy must actually be estimated/used, so flip it back on
    if ("rotate" %in% names(spcov_initial$initial)) {
      is_rotate_zero <- (!is.na(spcov_initial$initial[["rotate"]])) && spcov_initial$initial[["rotate"]] == 0
      is_rotate_known <- spcov_initial$is_known[["rotate"]]
      if (!is_rotate_zero || !is_rotate_known) {
        anisotropy <- TRUE
      }
    }

    if ("scale" %in% names(spcov_initial$initial)) {
      is_scale_one <- (!is.na(spcov_initial$initial[["scale"]])) && spcov_initial$initial[["scale"]] == 1
      is_scale_known <- spcov_initial$is_known[["scale"]]
      if (!is_scale_one || !is_scale_known) {
        anisotropy <- TRUE
      }
    }
  }
  anisotropy
}
