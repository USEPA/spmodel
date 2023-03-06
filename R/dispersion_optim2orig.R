dispersion_optim2orig <- function(dispersion_orig2optim, par) {
  par <- par[length(par)] # dispersion is the last element
  fill_optim_par_val <- fill_optim_par(dispersion_orig2optim, par)
  fill_orig_val <- exp(fill_optim_par_val)
  names(fill_orig_val) <- "dispersion"
  fill_orig_val
}
