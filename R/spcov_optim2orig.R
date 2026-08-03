#' Transform spatial covariance parameters from optim to original scale
#'
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param par Parameters to optimize over
#' @param spcov_profiled Is spatial profiling used?
#'
#' @return Spatial covariance parameters on the original scale
#'
#' @noRd
spcov_optim2orig <- function(spcov_orig2optim, par, spcov_profiled, ...) {
  # inverse of spcov_orig2optim(): maps the unconstrained parameters used
  # internally by optim() back to the original, constrained covariance-
  # parameter scale (e.g., variances back from log-scale, proportions/angles
  # back from logit-scale). Dispatches by class, one method per covariance type
  UseMethod("spcov_optim2orig", spcov_orig2optim)
}

#' @export
spcov_optim2orig.exponential <- function(spcov_orig2optim, par, spcov_profiled, data_object, ...) { # data object not used for geostatistical models
  # combine the values actually being optimized (par) with any values fixed
  # as "known" into the full named parameter vector expected below
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  de_ie <- optim2orig_de_ie(fill_optim_par_val, spcov_profiled)
  range <- optim2orig_range(fill_optim_par_val, data_object)
  aniso <- optim2orig_anisotropy(fill_optim_par_val)

  fill_orig_val <- c(de_ie, range = range, aniso)
}

#' @export
# these covariance types share the exponential method because they differ
# only in their spatial correlation form R(h), not in how de/ie/range/
# rotate/scale are parameterized for optimization
spcov_optim2orig.spherical <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.gaussian <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.triangular <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.circular <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.none <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.ie <- spcov_optim2orig.none
#' @export
spcov_optim2orig.cubic <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.pentaspherical <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.cosine <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.wave <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.jbessel <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.gravity <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.rquad <- spcov_optim2orig.exponential
#' @export
spcov_optim2orig.magnetic <- spcov_optim2orig.exponential

#' @export
spcov_optim2orig.matern <- function(spcov_orig2optim, par, spcov_profiled, data_object, ...) {
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  de_ie <- optim2orig_de_ie(fill_optim_par_val, spcov_profiled)

  # extra (matern smoothness) is bounded to [1/5, 5]; expit maps the
  # unconstrained optimizer value to (0, 1), then it's rescaled into range
  extra_t <- expit(fill_optim_par_val[["extra_logodds"]])
  # fix to be in [1/5, 5]
  extra <- extra_t * (5 - 1 / 5) + 1 / 5

  range <- optim2orig_range(fill_optim_par_val, data_object)
  aniso <- optim2orig_anisotropy(fill_optim_par_val)

  fill_orig_val <- c(de_ie, range = range, extra = extra, aniso)
}

#' @export
spcov_optim2orig.cauchy <- function(spcov_orig2optim, par, spcov_profiled, data_object, ...) {
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  de_ie <- optim2orig_de_ie(fill_optim_par_val, spcov_profiled)

  extra <- exp(fill_optim_par_val[["extra_log"]])

  range <- optim2orig_range(fill_optim_par_val, data_object)
  aniso <- optim2orig_anisotropy(fill_optim_par_val)

  fill_orig_val <- c(de_ie, range = range, extra = extra, aniso)
}

#' @export
spcov_optim2orig.pexponential <- function(spcov_orig2optim, par, spcov_profiled, data_object, ...) {
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  de_ie <- optim2orig_de_ie(fill_optim_par_val, spcov_profiled)

  # extra (pexponential power) is bounded to (0, 2]; expit gives (0, 1), then
  # doubled to land in the valid range
  extra_half <- expit(fill_optim_par_val[["extra_logodds"]])
  extra <- 2 * extra_half

  range <- optim2orig_range(fill_optim_par_val, data_object)
  aniso <- optim2orig_anisotropy(fill_optim_par_val)

  fill_orig_val <- c(de_ie, range = range, extra = extra, aniso)
}

#' @export
spcov_optim2orig.car <- function(spcov_orig2optim, par, spcov_profiled, data_object, ...) {
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  if (spcov_profiled) {
    # car profiling additionally splits off extra's share of the total
    # variance (in addition to ie's share), since car has a third variance
    # component (extra, for unconnected observations)
    ie_prop <- expit(fill_optim_par_val[["ie_prop_logodds"]])
    extra_prop <- expit(fill_optim_par_val[["extra_prop_logodds"]])
    de <- (1 - ie_prop) * (1 - extra_prop)
    ie <- ie_prop * (1 - extra_prop)
    extra <- extra_prop
  } else {
    de <- exp(fill_optim_par_val[["de_log"]])
    ie <- exp(fill_optim_par_val[["ie_log"]])
    extra <- exp(fill_optim_par_val[["extra_log"]])
  }

  # range (spatial autocorrelation) is valid only between the reciprocal
  # eigenvalue bounds of W (rho_lb, rho_ub); optimize on logit scale then map
  # (0, 1) into that valid interval
  range <- expit(fill_optim_par_val[["range_logodds"]])
  range <- range * (data_object$rho_ub - data_object$rho_lb) + data_object$rho_lb # scale to proper value

  fill_orig_val <- c(de = de, ie = ie, range = range, extra = extra)
}
#' @export
spcov_optim2orig.sar <- spcov_optim2orig.car
