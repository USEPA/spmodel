#' Set NA defaults for spcov_iniital objects
#'
#' @param spcov_initial A \code{spcov_initial object}
#'
#' @return A spcov_initial object with NA values when relevant (which are replaced later)
#' -- values are NA if they want us to pick those initial values
#'
#' @noRd
spcov_initial_NA <- function(spcov_initial, anisotropy = FALSE, is_W_connected = NULL) {
  # three parameter family
  if (inherits(spcov_initial, c("exponential", "spherical", "gaussian", "triangular", "circular", "cubic", "pentaspherical", "cosine", "wave", "jbessel", "gravity", "rquad", "magnetic"))) {
    spcov_names <- c("de", "ie", "range", "rotate", "scale")
    if (anisotropy) {
      spcov_val_default <- c(de = NA, ie = NA, range = NA, rotate = NA, scale = NA)
      spcov_known_default <- c(de = FALSE, ie = FALSE, range = FALSE, rotate = FALSE, scale = FALSE)
    } else {
      spcov_val_default <- c(de = NA, ie = NA, range = NA, rotate = 0, scale = 1)
      spcov_known_default <- c(de = FALSE, ie = FALSE, range = FALSE, rotate = TRUE, scale = TRUE)
    }
  } else if (inherits(spcov_initial, c("matern", "cauchy", "pexponential"))) { # 4 parameter family geo
    spcov_names <- c("de", "ie", "range", "extra", "rotate", "scale")
    if (anisotropy) {
      spcov_val_default <- c(de = NA, ie = NA, range = NA, extra = NA, rotate = NA, scale = NA)
      spcov_known_default <- c(de = FALSE, ie = FALSE, range = FALSE, extra = FALSE, rotate = FALSE, scale = FALSE)
    } else {
      spcov_val_default <- c(de = NA, ie = NA, range = NA, extra = NA, rotate = 0, scale = 1)
      spcov_known_default <- c(de = FALSE, ie = FALSE, range = FALSE, extra = FALSE, rotate = TRUE, scale = TRUE)
    }
  } else if (inherits(spcov_initial, c("car", "sar"))) { # 4 parameter family ar
    if (is_W_connected) {
      spcov_names <- c("de", "ie", "range", "extra")
      spcov_val_default <- c(de = NA, ie = 0, range = NA, extra = 0)
      spcov_known_default <- c(de = FALSE, ie = TRUE, range = FALSE, extra = TRUE)
    } else {
      spcov_names <- c("de", "ie", "range", "extra")
      spcov_val_default <- c(de = NA, ie = 0, range = NA, extra = NA)
      spcov_known_default <- c(de = FALSE, ie = TRUE, range = FALSE, extra = FALSE)
    }
  } else if (inherits(spcov_initial, c("none", "ie"))) {
    spcov_names <- c("de", "ie", "range", "rotate", "scale")
    spcov_val_default <- c(de = 0, ie = NA, range = Inf, rotate = 0, scale = 1)
    spcov_known_default <- c(de = TRUE, ie = FALSE, range = TRUE, rotate = TRUE, scale = TRUE)
  }

  # find names not in initial
  spcov_out <- setdiff(spcov_names, names(spcov_initial$initial))

  # put in values not in initial
  spcov_initial$initial[spcov_out] <- spcov_val_default[spcov_out]
  # put in is_known not in initial
  spcov_initial$is_known[spcov_out] <- spcov_known_default[spcov_out]
  # reorder names
  if (inherits(spcov_initial, c("none", "ie"))) {
    # reset if none covariance
    spcov_initial$initial[c("de", "range", "rotate", "scale")] <- spcov_val_default[c("de", "range", "rotate", "scale")]
    # put in is_known not in initial
    spcov_initial$is_known[c("de", "range", "rotate", "scale")] <- spcov_known_default[c("de", "range", "rotate", "scale")]
  }
  spcov_initial$initial <- spcov_initial$initial[spcov_names]
  spcov_initial$is_known <- spcov_initial$is_known[spcov_names]
  # return spcov_initial
  spcov_initial
}
