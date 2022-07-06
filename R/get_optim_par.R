#' Get parameters to optimize over in optim (remove known parameters)
#'
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param randcov_orig2optim A \code{randcov_orig2optim} object
#'
#' @return The parameters to optimize over in optim
#'
#' @noRd
get_optim_par <- function(spcov_orig2optim, randcov_orig2optim = NULL) {
  if (is.null(randcov_orig2optim)) {
    par <- spcov_orig2optim$value[!spcov_orig2optim$is_known]
  } else {
    spcov_pars <- spcov_orig2optim$value[!spcov_orig2optim$is_known]
    randcov_pars <- randcov_orig2optim$value[!randcov_orig2optim$is_known]
    par <- c(spcov_pars, randcov_pars)
  }
  par
}
