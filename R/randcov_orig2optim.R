#' Transform random effects from original to optim scale
#'
#' @param randcov_initial A \code{randcov_initial} object
#'
#' @return A \code{randcov_initial} list on the optim scale
#'
#' @noRd
randcov_orig2optim <- function(randcov_initial, randcov_profiled = NULL, spcov_initial = NULL) {
  if (is.null(randcov_initial)) {
    randcov_initial_list_val <- NULL
  } else {
    if (!is.null(randcov_profiled) && randcov_profiled) {
      # only done if all random effects unknown
      s2r <- randcov_initial$initial
      s2 <- sum(spcov_initial$initial[["de"]], spcov_initial$initial[["ie"]], s2r)
      v_r <- sum(s2r) / s2
      randcov_orig2optim_val <- c(v_r = v_r)
      n_randcov <- length(s2r)
      if (n_randcov > 1) {
        for (i in seq(2, n_randcov)) {
          randcov_orig2optim_val[i] <- sum(s2r[seq(i, n_randcov)]) / sum(s2r[seq(i - 1, n_randcov)])
          names(randcov_orig2optim_val)[i] <- paste("v", i, sep = "_")
        }
      }
      randcov_orig2optim_val <- log(randcov_orig2optim_val / (1 - randcov_orig2optim_val))
      randcov_orig2optim_is_known <- rep(FALSE, length(randcov_orig2optim_val))
      names(randcov_orig2optim_is_known) <- names(randcov_initial$is_known) # non profiled to keep for later
      randcov_initial_list_val <- list(
        value = randcov_orig2optim_val,
        is_known = randcov_orig2optim_is_known,
        n_est = sum(!randcov_orig2optim_is_known)
      )
    } else {
      randcov_orig2optim_val <- log(randcov_initial$initial)
      names(randcov_orig2optim_val) <- paste(names(randcov_initial$initial), "log", sep = "_")
      randcov_orig2optim_is_known <- randcov_initial$is_known
      names(randcov_orig2optim_is_known) <- paste(names(randcov_initial$is_known), "log", sep = "_")
      randcov_initial_list_val <- list(
        value = randcov_orig2optim_val,
        is_known = randcov_orig2optim_is_known,
        n_est = sum(!randcov_orig2optim_is_known)
      )
    }
  }
  randcov_initial_list_val
}
