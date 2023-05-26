dispersion_initial_NA <- function(dispersion_initial, data_object) {
  dispersion_names <- c("dispersion")
  if (data_object$family %in% c("poisson", "binomial")) {
    dispersion_initial <- dispersion_initial(data_object$family, 1, known = "dispersion")
  } else {
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
