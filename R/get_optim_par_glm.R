get_optim_par_glm <- function(spcov_orig2optim, dispersion_orig2optim, randcov_orig2optim = NULL) {
  if (is.null(randcov_orig2optim)) {
    par <- spcov_orig2optim$value[!spcov_orig2optim$is_known]
  } else {
    spcov_pars <- spcov_orig2optim$value[!spcov_orig2optim$is_known]
    randcov_pars <- randcov_orig2optim$value[!randcov_orig2optim$is_known]
    par <- c(spcov_pars, randcov_pars)
  }
  dispersion_pars <- dispersion_orig2optim$value[!dispersion_orig2optim$is_known]
  par <- c(par, dispersion_pars)
  par
}
