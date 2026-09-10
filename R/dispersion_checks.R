#' Validate a user-supplied dispersion value against the response family
#'
#' @param family The response family
#' @param dispersion A dispersion value (or \code{NULL})
#'
#' @return Nothing; errors if \code{dispersion} is fixed at a value other than
#'   one for the binomial or Poisson families (whose dispersion must be one)
#'
#' @noRd
dispersion_checks <- function(family, dispersion) {
  # dispersion can't be missing
  # binomial and poisson have a variance function fully determined by the mean,
  # so their dispersion parameter is not identifiable and must stay fixed at one
  if (!is.null(dispersion) && dispersion != 1 && family %in% c("binomial", "poisson")) {
    stop(paste(family, "dispersion parameter must be fixed at one."), call. = FALSE)
  }
}
