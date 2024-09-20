spcov_initial_NA_glm <- function(family, spcov_initial, anisotropy = FALSE, is_W_connected = NULL) {
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy, is_W_connected)
  spcov_initial_NA_glm_val <- spcov_initial_NA_val

  if (inherits(spcov_initial_NA_glm_val, "none")) {
    spcov_initial_NA_glm_val$initial["ie"] <- 0
    spcov_initial_NA_glm_val$is_known["ie"] <- TRUE
  }

  # fix ie can be confounded with dispersion (problems with ie error)
  # if (family %in% c("nbinomial", "Gamma", "inverse.gaussian", "beta")) {
  #   if (is.na(spcov_initial_NA_glm_val$initial[["ie"]])) {
  #     if (inherits(spcov_initial, "none")) {
  #       spcov_initial_NA_glm_val$initial[["ie"]] <- 1e-4
  #     } else {
  #       spcov_initial_NA_glm_val$initial[["ie"]] <- 1e-4 # or 0
  #     }
  #     spcov_initial_NA_glm_val$is_known[["ie"]] <- TRUE
  #   }
  # }
  spcov_initial_NA_glm_val
}
