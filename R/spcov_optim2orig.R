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
  UseMethod("spcov_optim2orig", spcov_orig2optim)
}

#' @export
spcov_optim2orig.exponential <- function(spcov_orig2optim, par, spcov_profiled, data_object) { # data object not used for geostatistical models
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  if (spcov_profiled) {
    ie_prop <- expit(fill_optim_par_val[["ie_prop_logodds"]])
    de <- 1 - ie_prop
    ie <- ie_prop
  } else {
    de <- exp(fill_optim_par_val[["de_log"]])
    ie <- exp(fill_optim_par_val[["ie_log"]])
  }

  range <- exp(fill_optim_par_val[["range_log"]])
  rotate <- pi * expit(fill_optim_par_val[["rotate_logodds"]])
  scale <- expit(fill_optim_par_val[["scale_logodds"]])


  fill_orig_val <- c(de = de, ie = ie, range = range, rotate = rotate, scale = scale)
}

#' @export
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
spcov_optim2orig.matern <- function(spcov_orig2optim, par, spcov_profiled, data_object) {
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  if (spcov_profiled) {
    ie_prop <- expit(fill_optim_par_val[["ie_prop_logodds"]])
    de <- 1 - ie_prop
    ie <- ie_prop
  } else {
    de <- exp(fill_optim_par_val[["de_log"]])
    ie <- exp(fill_optim_par_val[["ie_log"]])
  }


  extra_t <- expit(fill_optim_par_val[["extra_logodds"]])
  # fix to be in [1/5, 5]
  extra <- extra_t * (5 - 1 / 5) + 1 / 5
  range <- exp(fill_optim_par_val[["range_log"]])
  rotate <- pi * expit(fill_optim_par_val[["rotate_logodds"]])
  scale <- expit(fill_optim_par_val[["scale_logodds"]])


  fill_orig_val <- c(de = de, ie = ie, range = range, extra = extra, rotate = rotate, scale = scale)
}

#' @export
spcov_optim2orig.cauchy <- function(spcov_orig2optim, par, spcov_profiled, data_object) {
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  if (spcov_profiled) {
    ie_prop <- expit(fill_optim_par_val[["ie_prop_logodds"]])
    de <- 1 - ie_prop
    ie <- ie_prop
  } else {
    de <- exp(fill_optim_par_val[["de_log"]])
    ie <- exp(fill_optim_par_val[["ie_log"]])
  }


  extra <- exp(fill_optim_par_val[["extra_log"]])
  range <- exp(fill_optim_par_val[["range_log"]])
  rotate <- pi * expit(fill_optim_par_val[["rotate_logodds"]])
  scale <- expit(fill_optim_par_val[["scale_logodds"]])


  fill_orig_val <- c(de = de, ie = ie, range = range, extra = extra, rotate = rotate, scale = scale)
}

#' @export
spcov_optim2orig.pexponential <- function(spcov_orig2optim, par, spcov_profiled, data_object) {
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  if (spcov_profiled) {
    ie_prop <- expit(fill_optim_par_val[["ie_prop_logodds"]])
    de <- 1 - ie_prop
    ie <- ie_prop
  } else {
    de <- exp(fill_optim_par_val[["de_log"]])
    ie <- exp(fill_optim_par_val[["ie_log"]])
  }

  extra_half <- expit(fill_optim_par_val[["extra_logodds"]])
  extra <- 2 * extra_half
  range <- exp(fill_optim_par_val[["range_log"]])
  rotate <- pi * expit(fill_optim_par_val[["rotate_logodds"]])
  scale <- expit(fill_optim_par_val[["scale_logodds"]])


  fill_orig_val <- c(de = de, ie = ie, range = range, extra = extra, rotate = rotate, scale = scale)
}

#' @export
spcov_optim2orig.car <- function(spcov_orig2optim, par, spcov_profiled, data_object) {
  fill_optim_par_val <- fill_optim_par(spcov_orig2optim, par[seq(1, spcov_orig2optim$n_est)])

  if (spcov_profiled) {
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

  range <- expit(fill_optim_par_val[["range_logodds"]])
  range <- range * (data_object$rho_ub - data_object$rho_lb) + data_object$rho_lb # scale to proper value

  fill_orig_val <- c(de = de, ie = ie, range = range, extra = extra)
}
#' @export
spcov_optim2orig.sar <- spcov_optim2orig.car
