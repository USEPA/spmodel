#' Fill random effect parameters with NA's and known FALSE if specified in formula
#'   but notgiven as iniital value
#'
#' @param randcov_initial A \code{randcov_initial} object
#' @param randcov_names Names of random effects
#'
#' @return A \code{randcov_initial} object with appropriate NA's
#'
#' @noRd
randcov_initial_NA <- function(randcov_initial = NULL, randcov_names) {
  if (is.null(randcov_initial)) {
    randcov_initial <- NULL
  } else {
    randcov_names <- randcov_names
    randcov_val_default <- rep(NA, length = length(randcov_names))
    names(randcov_val_default) <- randcov_names
    randcov_known_default <- rep(FALSE, length = length(randcov_names))
    names(randcov_known_default) <- randcov_names
    # find names not in initial
    randcov_out <- setdiff(randcov_names, names(randcov_initial$initial))
    # put in values not in initial
    randcov_initial$initial[randcov_out] <- randcov_val_default[randcov_out]
    # put in is_known not in initial
    randcov_initial$is_known[randcov_out] <- randcov_known_default[randcov_out]
    # reorder names
    randcov_initial$initial <- randcov_initial$initial[randcov_names]
    randcov_initial$is_known <- randcov_initial$is_known[randcov_names]
  }

  # return randcov_initial
  randcov_initial
}
