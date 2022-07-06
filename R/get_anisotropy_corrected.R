get_anisotropy_corrected <- function(anisotropy, spcov_initial) {
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
