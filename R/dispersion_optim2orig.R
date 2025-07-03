dispersion_optim2orig <- function(dispersion_orig2optim, par) {
  dispersion_is_known <- dispersion_orig2optim$is_known

  # don't actually use fill_optim_par function
  if (dispersion_is_known) {
    fill_optim_par_val <- dispersion_orig2optim$value
    fill_orig_val <- exp(fill_optim_par_val)
    names(fill_orig_val) <- "dispersion"
    new_par <- par # don't need to remove anything
  } else {
    fill_optim_par_val <- par[length(par)] # dispersion is the last element
    fill_orig_val <- exp(fill_optim_par_val)
    # cap lower and upper dispersion values for numeric stability
    fill_orig_val <- pmax(1e-8, fill_orig_val)
    fill_orig_val <- pmin(1e8, fill_orig_val)
    names(fill_orig_val) <- "dispersion"
    new_par <- par[-length(par)] # remove dispersion parameter
  }
  list(fill_orig_val = fill_orig_val, new_par = new_par)
}
