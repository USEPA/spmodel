dispersion_initial <- function(family, dispersion, known) {

  # fix family
  if (is.symbol(substitute(family))) { # or is.language
    family <- deparse1(substitute(family))
  }

  # set defaults
  if (missing(dispersion)) dispersion <- NULL

  dispersion_checks(family, dispersion)

  dispersion_params_given <- c(dispersion = unname(dispersion))

  if (missing(known)) {
    is_known <- rep(FALSE, length(dispersion_params_given))
  } else {
    if (identical(known, "given")) {
      is_known <- rep(TRUE, length(dispersion_params_given))
    } else {
      is_known <- names(dispersion_params_given) %in% known
    }
  }
  names(is_known) <- names(dispersion_params_given)

  # error if NA and known
  dispersion_NA <- which(is.na(dispersion_params_given))
  if (any(is_known[dispersion_NA])) {
    stop("dispersion_initial values cannot be NA and known.", call. = FALSE)
  }

  new_dispersion_initial <- structure(list(initial = dispersion_params_given, is_known = is_known),
                                      class = family)
  new_dispersion_initial

}
