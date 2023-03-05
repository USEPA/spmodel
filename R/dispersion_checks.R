dispersion_checks <- function(family, dispersion) {

  # family must be a character here
  family_valid <- c("binomial", "poisson", "nbinomial", "Gamma")
  if (!(family %in% family_valid)) {
    stop(paste(family, " is not a valid glm family."), call. = FALSE)
  }
  # dispersion can't be missing
  if (!is.null(dispersion) && dispersion != 1 && family %in% c("binomial", "poisson")) {
    stop(paste(family, "dispersion parameter must be fixed at one."), call. = FALSE)
  }
}
