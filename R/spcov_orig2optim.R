#' Transform spatial covariance parameters from original to optim scale
#'
#' @param spcov_initial An \code{spcov_initial} object
#' @param spcov_profiled Is spatial profiling used?
#'
#' @return Covariance parameters on the optimi scale
#'
#' @noRd
spcov_orig2optim <- function(spcov_initial, spcov_profiled, ...) {
  # inverse of spcov_optim2orig(): maps covariance parameters on their
  # natural (constrained) scale to an unconstrained scale (log for positive
  # variances, logit for proportions/bounded parameters) so an unconstrained
  # optimizer like optim() can search freely without violating parameter
  # constraints. Dispatches by class, one method per covariance type
  UseMethod("spcov_orig2optim", spcov_initial)
}

#' @export
spcov_orig2optim.exponential <- function(spcov_initial, spcov_profiled, data_object, ...) { # data object not used for geostatistical models
  # exponential's ie_prop_logodds is_known reflects real conditional logic
  # (TRUE if de/ie both known, or ie known-and-zero) -- see orig2optim_de_ie()
  de_ie <- orig2optim_de_ie(spcov_initial, spcov_profiled, smart_is_known = TRUE)
  rng <- orig2optim_range(spcov_initial, data_object)
  aniso <- orig2optim_anisotropy(spcov_initial)

  value <- c(de_ie$value, rng$value, aniso$value)
  is_known <- c(de_ie$is_known, rng$is_known, aniso$is_known)

  # class = class(spcov_initial) lets spcov_optim2orig() dispatch back to the
  # matching inverse-transform method later
  finalize_orig2optim(value, is_known, spcov_initial)
}

#' @export
spcov_orig2optim.spherical <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.gaussian <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.triangular <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.circular <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.none <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.ie <- spcov_orig2optim.none
#' @export
spcov_orig2optim.cubic <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.pentaspherical <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.cosine <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.wave <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.jbessel <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.gravity <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.rquad <- spcov_orig2optim.exponential
#' @export
spcov_orig2optim.magnetic <- spcov_orig2optim.exponential

#' @export
spcov_orig2optim.matern <- function(spcov_initial, spcov_profiled, data_object, ...) {
  # matern's ie_prop_logodds is_known is always FALSE, unlike exponential's --
  # see orig2optim_de_ie()
  de_ie <- orig2optim_de_ie(spcov_initial, spcov_profiled, smart_is_known = FALSE)
  rng <- orig2optim_range(spcov_initial, data_object)

  # fix in [1/5, 5]
  # matern smoothness (extra) is bounded to [1/5, 5]; rescale to [0, 1] first
  # so logit gives a valid unconstrained value
  extra <- (spcov_initial$initial[["extra"]] - 1 / 5) / (5 - 1 / 5) # to be in [0, 1]
  extra_logodds <- logit(extra)

  aniso <- orig2optim_anisotropy(spcov_initial)

  value <- c(de_ie$value, rng$value, extra_logodds = extra_logodds, aniso$value)
  is_known <- c(de_ie$is_known, rng$is_known, extra_logodds = spcov_initial$is_known[["extra"]], aniso$is_known)

  finalize_orig2optim(value, is_known, spcov_initial)
}
#' @export
spcov_orig2optim.cauchy <- function(spcov_initial, spcov_profiled, data_object, ...) {
  # cauchy's ie_prop_logodds is_known is always FALSE, unlike exponential's --
  # see orig2optim_de_ie()
  de_ie <- orig2optim_de_ie(spcov_initial, spcov_profiled, smart_is_known = FALSE)
  rng <- orig2optim_range(spcov_initial, data_object)

  # extra p log
  extra <- spcov_initial$initial[["extra"]]
  extra_log <- log(extra)

  aniso <- orig2optim_anisotropy(spcov_initial)

  value <- c(de_ie$value, rng$value, extra_log = extra_log, aniso$value)
  is_known <- c(de_ie$is_known, rng$is_known, extra_log = spcov_initial$is_known[["extra"]], aniso$is_known)

  finalize_orig2optim(value, is_known, spcov_initial)
}
#' @export
spcov_orig2optim.pexponential <- function(spcov_initial, spcov_profiled, data_object, ...) {
  # pexponential's ie_prop_logodds is_known is always FALSE, unlike
  # exponential's -- see orig2optim_de_ie()
  de_ie <- orig2optim_de_ie(spcov_initial, spcov_profiled, smart_is_known = FALSE)
  rng <- orig2optim_range(spcov_initial, data_object)

  # extra p logodds (for now)
  # pexponential's extra is bounded to (0, 2]; halve it to [0, 1] before logit
  extra <- spcov_initial$initial[["extra"]]
  extra_half <- extra / 2 # because maximum value is 2
  extra_logodds <- logit(extra_half)

  aniso <- orig2optim_anisotropy(spcov_initial)

  value <- c(de_ie$value, rng$value, extra_logodds = extra_logodds, aniso$value)
  is_known <- c(de_ie$is_known, rng$is_known, extra_logodds = spcov_initial$is_known[["extra"]], aniso$is_known)

  finalize_orig2optim(value, is_known, spcov_initial)
}

#' @export
spcov_orig2optim.car <- function(spcov_initial, spcov_profiled, data_object, ...) {
  # are variance parameters spcov_profiled
  if (spcov_profiled) {
    # profiling is only supported when extra is fixed at 0 (i.e., no
    # unconnected observations) -- with a nonzero/estimated extra there is no
    # single total-variance term to profile out, hence the stop() below
    if (spcov_initial$initial[["extra"]] == 0 && spcov_initial$is_known[["extra"]]) {
      # log odds
      ie_prop <- spcov_initial$initial[["ie"]] / (spcov_initial$initial[["de"]] + spcov_initial$initial[["ie"]])
      ie_prop_logodds <- logit(ie_prop)
      if (spcov_initial$is_known[["de"]] && spcov_initial$is_known[["ie"]]) {
        ie_prop_logodds_is_known <- TRUE
        # } else if (spcov_initial$is_known[["de"]] && spcov_initial$initial[["de"]] == 0) {
        #   ie_prop_logodds_is_known <- TRUE # not needed here because iid would be called
      } else if (spcov_initial$is_known[["ie"]] && spcov_initial$initial[["ie"]] == 0) {
        ie_prop_logodds_is_known <- TRUE
      } else {
        ie_prop_logodds_is_known <- FALSE
      }

      extra_prop <- 0
      extra_prop_logodds <- logit(extra_prop)
      extra_prop_logodds_is_known <- TRUE

      spcov_orig2optim_val <- c(ie_prop_logodds = ie_prop_logodds, extra_prop_logodds = extra_prop_logodds)
      spcov_orig2optim_is_known <- c(ie_prop_logodds = ie_prop_logodds_is_known, extra_prop_logodds = extra_prop_logodds_is_known)
    } else {
      stop("Currently, profiling and extra cannot be used simultaneously")
    }
  } else { # log
    de <- spcov_initial$initial[["de"]]
    de_log <- log(spcov_initial$initial[["de"]])
    ie <- spcov_initial$initial[["ie"]]
    ie_log <- log(spcov_initial$initial[["ie"]])
    extra <- spcov_initial$initial[["extra"]]
    extra_log <- log(spcov_initial$initial[["extra"]])
    spcov_orig2optim_val <- c(de_log = de_log, ie_log = ie_log, extra_log = extra_log)
    spcov_orig2optim_is_known <- c(
      de_log = spcov_initial$is_known[["de"]],
      ie_log = spcov_initial$is_known[["ie"]],
      extra_log = spcov_initial$is_known[["extra"]]
    )
  }

  # range changes based on type
  # car/sar range (spatial autocorrelation) is only valid between the
  # reciprocal eigenvalue bounds of W (rho_lb, rho_ub); rescale to [0, 1]
  # before the logit transform
  range <- spcov_initial$initial[["range"]]
  range <- (range - data_object$rho_lb) / (data_object$rho_ub - data_object$rho_lb) # scale to 0,1
  # range <- (range + 1) / 2 # (from -1, 1 to 0, 2 to 0, 1)
  range_logodds <- logit(range)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, range_logodds = range_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_logodds = spcov_initial$is_known[["range"]])

  finalize_orig2optim(spcov_orig2optim_val, spcov_orig2optim_is_known, spcov_initial)
}
#' @export
spcov_orig2optim.sar <- spcov_orig2optim.car
