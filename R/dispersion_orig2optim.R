dispersion_orig2optim <- function(dispersion_initial) {
  dispersion_orig2optim_val <- log(dispersion_initial$initial)
  names(dispersion_orig2optim_val) <- "dispersion_log"
  dispersion_orig2optim_is_known <- dispersion_initial$is_known
  names(dispersion_orig2optim_is_known) <- "dispersion_log"
  dispersion_initial_list_val <- list(
    value = dispersion_orig2optim_val,
    is_known = dispersion_orig2optim_is_known
  )
  dispersion_initial_list_val
}
