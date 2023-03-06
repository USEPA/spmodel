dispersion_params <- function(family, dispersion) {

  # fix family
  if (is.symbol(substitute(family))) { # or is.language
    family <- deparse1(substitute(family))
  }

  dispersion_checks(family, dispersion)

  object <- c(dispersion = unname(dispersion))
  new_object <- structure(object, class = family)
  new_object
}
