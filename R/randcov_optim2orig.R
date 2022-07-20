#' Transform random effects from optim to original scale (noting length of spcov params scale)
#'
#' @param randcov_orig2optim A \code{randcov_orig2optim} object
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param par Optim parameters
#'
#' @return Paramters on the origin scale
#'
#' @noRd
randcov_optim2orig <- function(randcov_orig2optim = NULL, spcov_orig2optim,
                               par, randcov_profiled = NULL, spcov_optim2orig = NULL) {
  if (is.null(randcov_orig2optim)) {
    fill_orig_val <- NULL
  } else {
    par_NA <- par[seq(spcov_orig2optim$n_est + 1, length(par))] # if all random effects known
    par <- par_NA[!is.na(par_NA)]
    fill_optim_par_val <- fill_optim_par(randcov_orig2optim, par)
    if (!is.null(randcov_profiled) && randcov_profiled) {
      fill_orig_val <- 1 / (1 + (1 / exp(fill_optim_par_val)))
      v_r <- fill_orig_val[1]
      spcov_optim2orig[["de"]] <- (1 - v_r) * spcov_optim2orig[["de"]]
      spcov_optim2orig[["ie"]] <- (1 - v_r) * spcov_optim2orig[["ie"]]
      n_randcov <- length(fill_orig_val)
      if (n_randcov > 1) {
        if (n_randcov > 2) {
          fill_orig_val[n_randcov] <- prod(fill_orig_val[-1])
          for (i in seq(n_randcov - 1, 2)) {
            back_index <- seq(i + 1, n_randcov)
            fill_orig_val[i] <- prod(fill_orig_val[-c(1, back_index)]) - sum(fill_orig_val[back_index])
          }
        }
        fill_orig_val[1] <- 1 - sum(fill_orig_val[-1])
        fill_orig_val <- v_r * fill_orig_val
      }
      names(fill_orig_val) <- names(randcov_orig2optim$is_known)
      fill_orig_val <- list(fill_orig_val = fill_orig_val, spcov_optim2orig = spcov_optim2orig)
    } else {
      fill_orig_val <- exp(fill_optim_par_val)
      names(fill_orig_val) <- gsub("_log", "", names(randcov_orig2optim$value))
    }
  }
  fill_orig_val
}
