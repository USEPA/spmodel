# shared building blocks for spcov_optim2orig.R / spcov_orig2optim.R's
# exponential/matern/cauchy/pexponential methods -- extracted because these
# four bodies duplicate the same de/ie, range, and rotate/scale transform
# logic verbatim, differing only in their own inline `extra` parameter
# handling. car/sar are NOT covered here: their de/ie/range transforms are
# structurally different (a third variance component, rho_lb/rho_ub-bounded
# range, no anisotropy) and stay fully inline in spcov_orig2optim.car/
# spcov_optim2orig.car.

#' Recover de/ie from their optim-scale representation
#'
#' @param fill_optim_par_val A named numeric vector produced by
#'   \code{fill_optim_par()}, containing either \code{ie_prop_logodds} (when
#'   \code{spcov_profiled}) or \code{de_log}/\code{ie_log}
#' @param spcov_profiled Is spatial profiling used?
#'
#' @return A named numeric vector \code{c(de =, ie =)} on the original scale
#'
#' @noRd
optim2orig_de_ie <- function(fill_optim_par_val, spcov_profiled) {
  if (spcov_profiled) {
    ie_prop <- expit(fill_optim_par_val[["ie_prop_logodds"]])
    de <- 1 - ie_prop
    ie <- ie_prop
  } else {
    de <- exp(fill_optim_par_val[["de_log"]])
    ie <- exp(fill_optim_par_val[["ie_log"]])
  }
  c(de = de, ie = ie)
}

#' Recover range from its optim-scale representation
#'
#' @param fill_optim_par_val A named numeric vector produced by
#'   \code{fill_optim_par()}, containing either \code{range_logodds} (when
#'   \code{data_object$range_constrain}) or \code{range_log}
#' @param data_object A data object with elements \code{range_constrain} and
#'   \code{range_constrain_value}
#'
#' @return \code{range} on the original scale
#'
#' @noRd
optim2orig_range <- function(fill_optim_par_val, data_object) {
  if (data_object$range_constrain) {
    range <- expit(fill_optim_par_val[["range_logodds"]]) * data_object$range_constrain_value
  } else {
    range <- exp(fill_optim_par_val[["range_log"]])
  }
  range
}

#' Recover rotate/scale (anisotropy) from their optim-scale representation
#'
#' @param fill_optim_par_val A named numeric vector produced by
#'   \code{fill_optim_par()}, containing \code{rotate_logodds}/\code{scale_logodds}
#'
#' @return A named numeric vector \code{c(rotate =, scale =)} on the original scale
#'
#' @noRd
optim2orig_anisotropy <- function(fill_optim_par_val) {
  rotate <- pi * expit(fill_optim_par_val[["rotate_logodds"]])
  scale <- expit(fill_optim_par_val[["scale_logodds"]])
  c(rotate = rotate, scale = scale)
}

#' Transform de/ie to their optim-scale representation
#'
#' @param spcov_initial An \code{spcov_initial} object
#' @param spcov_profiled Is spatial profiling used?
#' @param smart_is_known If \code{TRUE} (exponential and its aliases),
#'   \code{ie_prop_logodds}'s \code{is_known} is \code{TRUE} whenever de and ie
#'   are both known, or ie is known and zero -- otherwise \code{FALSE}. If
#'   \code{FALSE} (matern/cauchy/pexponential), \code{ie_prop_logodds}'s
#'   \code{is_known} is always \code{FALSE} regardless of what's known. This is
#'   a genuine pre-existing difference between covariance types, not something
#'   to unify -- callers must state which behavior they're preserving.
#'
#' @return A list with elements \code{value} and \code{is_known}
#'
#' @noRd
orig2optim_de_ie <- function(spcov_initial, spcov_profiled, smart_is_known) {
  if (spcov_profiled) {
    ie_prop <- spcov_initial$initial[["ie"]] / (spcov_initial$initial[["de"]] + spcov_initial$initial[["ie"]])
    ie_prop_logodds <- logit(ie_prop)
    value <- c(ie_prop_logodds = ie_prop_logodds)
    if (smart_is_known) {
      if (spcov_initial$is_known[["de"]] && spcov_initial$is_known[["ie"]]) {
        ie_prop_logodds_is_known <- TRUE
      } else if (spcov_initial$is_known[["ie"]] && spcov_initial$initial[["ie"]] == 0) {
        ie_prop_logodds_is_known <- TRUE
      } else {
        ie_prop_logodds_is_known <- FALSE
      }
    } else {
      ie_prop_logodds_is_known <- FALSE
    }
    is_known <- c(ie_prop_logodds = ie_prop_logodds_is_known)
  } else {
    de_log <- log(spcov_initial$initial[["de"]])
    ie_log <- log(spcov_initial$initial[["ie"]])
    value <- c(de_log = de_log, ie_log = ie_log)
    is_known <- c(
      de_log = spcov_initial$is_known[["de"]],
      ie_log = spcov_initial$is_known[["ie"]]
    )
  }
  list(value = value, is_known = is_known)
}

#' Transform range to its optim-scale representation
#'
#' @param spcov_initial An \code{spcov_initial} object
#' @param data_object A data object with elements \code{range_constrain} and
#'   \code{range_constrain_value}
#'
#' @return A list with elements \code{value} and \code{is_known}
#'
#' @noRd
orig2optim_range <- function(spcov_initial, data_object) {
  range <- spcov_initial$initial[["range"]]
  if (data_object$range_constrain) {
    range_prop <- range / data_object$range_constrain_value
    range_logodds <- logit(range_prop)
    value <- c(range_logodds = range_logodds)
    is_known <- c(range_logodds = spcov_initial$is_known[["range"]])
  } else {
    range_log <- log(range)
    value <- c(range_log = range_log)
    is_known <- c(range_log = spcov_initial$is_known[["range"]])
  }
  list(value = value, is_known = is_known)
}

#' Transform rotate/scale (anisotropy) to their optim-scale representation
#'
#' @param spcov_initial An \code{spcov_initial} object
#'
#' @return A list with elements \code{value} and \code{is_known}
#'
#' @noRd
orig2optim_anisotropy <- function(spcov_initial) {
  rotate <- spcov_initial$initial[["rotate"]]
  rotate_prop <- rotate / pi
  rotate_logodds <- logit(rotate_prop)

  scale <- spcov_initial$initial[["scale"]]
  scale_logodds <- logit(scale)

  value <- c(rotate_logodds = rotate_logodds, scale_logodds = scale_logodds)
  is_known <- c(
    rotate_logodds = spcov_initial$is_known[["rotate"]],
    scale_logodds = spcov_initial$is_known[["scale"]]
  )
  list(value = value, is_known = is_known)
}

#' Clamp optim-scale values and wrap into an \code{spcov_orig2optim} object
#'
#' @param value A named numeric vector of optim-scale values
#' @param is_known A named logical vector, same length/order as \code{value}
#' @param spcov_initial An \code{spcov_initial} object (its class is carried
#'   over so \code{spcov_optim2orig()} dispatches back to the matching method)
#'
#' @return A \code{spcov_orig2optim} object: \code{list(value =, is_known =,
#'   n_est =)}, classed to match \code{spcov_initial}
#'
#' @noRd
finalize_orig2optim <- function(value, is_known, spcov_initial) {
  # clamp extreme log/logit values -- without this, values transformed back
  # via exp()/expit() during optimization could overflow/underflow, and the
  # clamp is skipped for known (fixed) parameters since those aren't searched
  value <- ifelse(value > 50 & !is_known, 50, value)
  value <- ifelse(value < -50 & !is_known, -50, value)

  spcov_orig2optim_val <- list(
    value = value,
    is_known = is_known,
    n_est = sum(!is_known)
  )

  structure(spcov_orig2optim_val, class = class(spcov_initial))
}
