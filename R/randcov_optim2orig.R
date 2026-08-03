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
    # random effect parameters occupy the tail of par, after the spcov
    # parameters; some entries may be NA if a random effect variance is known
    # (and thus not part of the vector optim() is searching over)
    par_NA <- par[seq(spcov_orig2optim$n_est + 1, length(par))] # if all random effects known
    par <- par_NA[!is.na(par_NA)]
    # merge the estimated values back in with any known (fixed) values to get
    # the full-length parameter vector on the optim scale
    fill_optim_par_val <- fill_optim_par(randcov_orig2optim, par)
    if (!is.null(randcov_profiled) && randcov_profiled) {
      # inverse-logit turns the optim-scale value back into a (0, 1) proportion;
      # this is the reverse of the transform in randcov_orig2optim()
      fill_orig_val <- 1 / (1 + (1 / exp(fill_optim_par_val)))
      # v_r is the share of total variance (de + ie + all random effects)
      # attributable to random effects as a whole
      v_r <- fill_orig_val[1]
      # de/ie make up the remaining (1 - v_r) share of total variance
      spcov_optim2orig[["de"]] <- (1 - v_r) * spcov_optim2orig[["de"]]
      spcov_optim2orig[["ie"]] <- (1 - v_r) * spcov_optim2orig[["ie"]]
      n_randcov <- length(fill_orig_val)
      if (n_randcov > 1) {
        if (n_randcov > 2) {
          # unwind the ratios (each v_i was the share of the
          # variance remaining after allocating to effects 1..i-1) back into
          # per-effect proportions, working from the last effect backwards
          fill_orig_val[n_randcov] <- prod(fill_orig_val[-1])
          for (i in seq(n_randcov - 1, 2)) {
            back_index <- seq(i + 1, n_randcov)
            fill_orig_val[i] <- prod(fill_orig_val[-c(1, back_index)]) - sum(fill_orig_val[back_index])
          }
        }
        # the first random effect gets whatever proportion is left over, then
        # everything is rescaled by v_r so the proportions sum to v_r (the
        # random effects' total share) rather than to 1
        fill_orig_val[1] <- 1 - sum(fill_orig_val[-1])
        fill_orig_val <- v_r * fill_orig_val
      }
      names(fill_orig_val) <- names(randcov_orig2optim$is_known)
      # de/ie were rescaled above, so they must travel back together with the
      # random effect variances for the caller to reconstruct a consistent set
      fill_orig_val <- list(fill_orig_val = fill_orig_val, spcov_optim2orig = spcov_optim2orig)
    } else {
      # variances are optimized on the log scale to keep them positive;
      # exponentiate to return to the original (variance) scale
      fill_orig_val <- exp(fill_optim_par_val)
      names(fill_orig_val) <- gsub("_log", "", names(randcov_orig2optim$value))
    }
  }
  fill_orig_val
}
