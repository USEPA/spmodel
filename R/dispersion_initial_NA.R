#' Fill in default (\code{NA}, unknown) dispersion initial values
#'
#' @param dispersion_initial A \code{dispersion_initial} object
#' @param data_object The data object
#'
#' @return A \code{dispersion_initial} object with any missing initial values
#'   and \code{is_known} indicators filled in with defaults (fixed at one for
#'   the binomial and Poisson families, otherwise unknown)
#'
#' @noRd
dispersion_initial_NA <- function(dispersion_initial, data_object) {
  dispersion_names <- c("dispersion")
  # poisson and binomial dispersion is not identifiable, so it is always fixed at one
  # regardless of what (if anything) the user supplied
  if (data_object$family %in% c("poisson", "binomial")) {
    dispersion_initial <- dispersion_initial(data_object$family, 1, known = "dispersion")
  } else {
    # any dispersion value the user did not specify defaults to NA (to be
    # estimated) and unknown, rather than erroring
    dispersion_val_default <- c(dispersion = NA)
    dispersion_known_default <- c(dispersion = FALSE)
    names_replace <- setdiff(dispersion_names, names(dispersion_initial$initial))
    dispersion_initial$initial[names_replace] <- dispersion_val_default[names_replace]
    dispersion_initial$is_known[names_replace] <- dispersion_known_default[names_replace]

    # reorder names
    dispersion_initial$initial <- dispersion_initial$initial[dispersion_names]
    dispersion_initial$is_known <- dispersion_initial$is_known[dispersion_names]
  }
  dispersion_initial
}
