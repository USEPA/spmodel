#' Transform spatial covariance parameters from original to optim scale
#'
#' @param spcov_initial An \code{spcov_initial} object
#' @param spcov_profiled Is spatial profiling used?
#'
#' @return Covariance parameters on the optimi scale
#'
#' @noRd
spcov_orig2optim <- function(spcov_initial, spcov_profiled, ...) {
  UseMethod("spcov_orig2optim", spcov_initial)
}

#' @export
spcov_orig2optim.exponential <- function(spcov_initial, spcov_profiled, data_object, ...) { # data object not used for geostatistical models
  # are variance parameters spcov_profiled
  if (spcov_profiled) { # log odds
    ie_prop <- spcov_initial$initial[["ie"]] / (spcov_initial$initial[["de"]] + spcov_initial$initial[["ie"]])
    ie_prop_logodds <- logit(ie_prop)
    spcov_orig2optim_val <- c(ie_prop_logodds = ie_prop_logodds)
    if (spcov_initial$is_known[["de"]] && spcov_initial$is_known[["ie"]]) {
      ie_prop_logodds_is_known <- TRUE
      # } else if (spcov_initial$is_known[["de"]] && spcov_initial$initial[["de"]] == 0) {
      #   ie_prop_logodds_is_known <- TRUE # not needed here because iid would be called
    } else if (spcov_initial$is_known[["ie"]] && spcov_initial$initial[["ie"]] == 0) {
      ie_prop_logodds_is_known <- TRUE
    } else {
      ie_prop_logodds_is_known <- FALSE
    }
    spcov_orig2optim_is_known <- c(ie_prop_logodds = ie_prop_logodds_is_known)
  } else { # log
    de <- spcov_initial$initial[["de"]]
    de_log <- log(spcov_initial$initial[["de"]])
    ie <- spcov_initial$initial[["ie"]]
    ie_log <- log(spcov_initial$initial[["ie"]])
    spcov_orig2optim_val <- c(de_log = de_log, ie_log = ie_log)
    spcov_orig2optim_is_known <- c(
      de_log = spcov_initial$is_known[["de"]],
      ie_log = spcov_initial$is_known[["ie"]]
    )
  }

  # range changes based on type
  range <- spcov_initial$initial[["range"]]
  if (data_object$range_constrain) {
    range_prop <- range / data_object$range_constrain_value
    range_logodds <- logit(range_prop)
    spcov_orig2optim_val <- c(spcov_orig2optim_val, range_logodds = range_logodds)
    spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_logodds = spcov_initial$is_known[["range"]])
  } else {
    range_log <- log(range)
    spcov_orig2optim_val <- c(spcov_orig2optim_val, range_log = range_log)
    spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_log = spcov_initial$is_known[["range"]])
  }

  # anisotropy parameters
  ## rotate (between 0 and pi radians)
  rotate <- spcov_initial$initial[["rotate"]]
  rotate_prop <- rotate / pi # used to be pi / 2
  rotate_logodds <- logit(rotate_prop)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, rotate_logodds = rotate_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, rotate_logodds = spcov_initial$is_known[["rotate"]])

  ## scale (between 0 and 1)
  scale <- spcov_initial$initial[["scale"]]
  scale_logodds <- logit(scale)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, scale_logodds = scale_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, scale_logodds = spcov_initial$is_known[["scale"]])

  # return covariance parameter vector
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val > 50 & !spcov_orig2optim_is_known, 50, spcov_orig2optim_val)
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val < -50 & !spcov_orig2optim_is_known, -50, spcov_orig2optim_val)

  # return list
  spcov_orig2optim_val <- list(
    value = spcov_orig2optim_val,
    is_known = spcov_orig2optim_is_known,
    n_est = sum(!spcov_orig2optim_is_known)
  )

  # give class to vector
  new_spcov_orig2optim_val <- structure(spcov_orig2optim_val, class = class(spcov_initial))
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
  # are variance parameters spcov_profiled
  if (spcov_profiled) { # log odds
    ie_prop <- spcov_initial$initial[["ie"]] / (spcov_initial$initial[["de"]] + spcov_initial$initial[["ie"]])
    ie_prop_logodds <- logit(ie_prop)
    spcov_orig2optim_val <- c(ie_prop_logodds = ie_prop_logodds)
    ie_prop_logodds_is_known <- FALSE
    spcov_orig2optim_is_known <- c(ie_prop_logodds = ie_prop_logodds_is_known)
  } else { # log
    de <- spcov_initial$initial[["de"]]
    de_log <- log(spcov_initial$initial[["de"]])
    ie <- spcov_initial$initial[["ie"]]
    ie_log <- log(spcov_initial$initial[["ie"]])
    spcov_orig2optim_val <- c(de_log = de_log, ie_log = ie_log)
    spcov_orig2optim_is_known <- c(
      de_log = spcov_initial$is_known[["de"]],
      ie_log = spcov_initial$is_known[["ie"]]
    )
  }

  # range changes based on type
  range <- spcov_initial$initial[["range"]]
  if (data_object$range_constrain) {
    range_prop <- range / data_object$range_constrain_value
    range_logodds <- logit(range_prop)
    spcov_orig2optim_val <- c(spcov_orig2optim_val, range_logodds = range_logodds)
    spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_logodds = spcov_initial$is_known[["range"]])
  } else {
    range_log <- log(range)
    spcov_orig2optim_val <- c(spcov_orig2optim_val, range_log = range_log)
    spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_log = spcov_initial$is_known[["range"]])
  }
  # range <- spcov_initial$initial[["range"]]
  # range_log <- log(range)
  # spcov_orig2optim_val <- c(spcov_orig2optim_val, range_log = range_log)
  # spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_log = spcov_initial$is_known[["range"]])

  # # extra p log (for now)
  # extra <- spcov_initial$initial[["extra"]]
  # extra_log <- log(extra)
  # spcov_orig2optim_val <- c(spcov_orig2optim_val, extra_log = extra_log)
  # spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, extra_log = spcov_initial$is_known[["extra"]])

  # fix in [1/5, 5]
  extra <- (spcov_initial$initial[["extra"]] - 1 / 5) / (5 - 1 / 5) # to be in [0, 1]
  extra_logodds <- logit(extra)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, extra_logodds = extra_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, extra_logodds = spcov_initial$is_known[["extra"]])

  # anisotropy parameters
  ## rotate (between 0 and pi radians)
  rotate <- spcov_initial$initial[["rotate"]]
  rotate_prop <- rotate / pi
  rotate_logodds <- logit(rotate_prop)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, rotate_logodds = rotate_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, rotate_logodds = spcov_initial$is_known[["rotate"]])

  ## scale (between 0 and 1)
  scale <- spcov_initial$initial[["scale"]]
  scale_logodds <- logit(scale)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, scale_logodds = scale_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, scale_logodds = spcov_initial$is_known[["scale"]])

  # return covariance parameter vector
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val > 50 & !spcov_orig2optim_is_known, 50, spcov_orig2optim_val)
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val < -50 & !spcov_orig2optim_is_known, -50, spcov_orig2optim_val)

  # return list
  spcov_orig2optim_val <- list(
    value = spcov_orig2optim_val,
    is_known = spcov_orig2optim_is_known,
    n_est = sum(!spcov_orig2optim_is_known)
  )

  # give class to vector
  new_spcov_orig2optim_val <- structure(spcov_orig2optim_val, class = class(spcov_initial))
}
#' @export
spcov_orig2optim.cauchy <- function(spcov_initial, spcov_profiled, data_object, ...) {
  # are variance parameters spcov_profiled
  if (spcov_profiled) { # log odds
    ie_prop <- spcov_initial$initial[["ie"]] / (spcov_initial$initial[["de"]] + spcov_initial$initial[["ie"]])
    ie_prop_logodds <- logit(ie_prop)
    spcov_orig2optim_val <- c(ie_prop_logodds = ie_prop_logodds)
    ie_prop_logodds_is_known <- FALSE
    spcov_orig2optim_is_known <- c(ie_prop_logodds = ie_prop_logodds_is_known)
  } else { # log
    de <- spcov_initial$initial[["de"]]
    de_log <- log(spcov_initial$initial[["de"]])
    ie <- spcov_initial$initial[["ie"]]
    ie_log <- log(spcov_initial$initial[["ie"]])
    spcov_orig2optim_val <- c(de_log = de_log, ie_log = ie_log)
    spcov_orig2optim_is_known <- c(
      de_log = spcov_initial$is_known[["de"]],
      ie_log = spcov_initial$is_known[["ie"]]
    )
  }

  # range changes based on type
  range <- spcov_initial$initial[["range"]]
  if (data_object$range_constrain) {
    range_prop <- range / data_object$range_constrain_value
    range_logodds <- logit(range_prop)
    spcov_orig2optim_val <- c(spcov_orig2optim_val, range_logodds = range_logodds)
    spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_logodds = spcov_initial$is_known[["range"]])
  } else {
    range_log <- log(range)
    spcov_orig2optim_val <- c(spcov_orig2optim_val, range_log = range_log)
    spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_log = spcov_initial$is_known[["range"]])
  }
  # range <- spcov_initial$initial[["range"]]
  # range_log <- log(range)
  # spcov_orig2optim_val <- c(spcov_orig2optim_val, range_log = range_log)
  # spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_log = spcov_initial$is_known[["range"]])

  # extra p log
  extra <- spcov_initial$initial[["extra"]]
  extra_log <- log(extra)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, extra_log = extra_log)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, extra_log = spcov_initial$is_known[["extra"]])

  # anisotropy parameters
  ## rotate (between 0 and pi radians)
  rotate <- spcov_initial$initial[["rotate"]]
  rotate_prop <- rotate / pi
  rotate_logodds <- logit(rotate_prop)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, rotate_logodds = rotate_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, rotate_logodds = spcov_initial$is_known[["rotate"]])

  ## scale (between 0 and 1)
  scale <- spcov_initial$initial[["scale"]]
  scale_logodds <- logit(scale)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, scale_logodds = scale_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, scale_logodds = spcov_initial$is_known[["scale"]])

  # return covariance parameter vector
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val > 50 & !spcov_orig2optim_is_known, 50, spcov_orig2optim_val)
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val < -50 & !spcov_orig2optim_is_known, -50, spcov_orig2optim_val)

  # return list
  spcov_orig2optim_val <- list(
    value = spcov_orig2optim_val,
    is_known = spcov_orig2optim_is_known,
    n_est = sum(!spcov_orig2optim_is_known)
  )

  # give class to vector
  new_spcov_orig2optim_val <- structure(spcov_orig2optim_val, class = class(spcov_initial))
}
#' @export
spcov_orig2optim.pexponential <- function(spcov_initial, spcov_profiled, data_object, ...) {
  # are variance parameters spcov_profiled
  if (spcov_profiled) { # log odds
    ie_prop <- spcov_initial$initial[["ie"]] / (spcov_initial$initial[["de"]] + spcov_initial$initial[["ie"]])
    ie_prop_logodds <- logit(ie_prop)
    spcov_orig2optim_val <- c(ie_prop_logodds = ie_prop_logodds)
    ie_prop_logodds_is_known <- FALSE
    spcov_orig2optim_is_known <- c(ie_prop_logodds = ie_prop_logodds_is_known)
  } else { # log
    de <- spcov_initial$initial[["de"]]
    de_log <- log(spcov_initial$initial[["de"]])
    ie <- spcov_initial$initial[["ie"]]
    ie_log <- log(spcov_initial$initial[["ie"]])
    spcov_orig2optim_val <- c(de_log = de_log, ie_log = ie_log)
    spcov_orig2optim_is_known <- c(
      de_log = spcov_initial$is_known[["de"]],
      ie_log = spcov_initial$is_known[["ie"]]
    )
  }

  # range changes based on type
  range <- spcov_initial$initial[["range"]]
  if (data_object$range_constrain) {
    range_prop <- range / data_object$range_constrain_value
    range_logodds <- logit(range_prop)
    spcov_orig2optim_val <- c(spcov_orig2optim_val, range_logodds = range_logodds)
    spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_logodds = spcov_initial$is_known[["range"]])
  } else {
    range_log <- log(range)
    spcov_orig2optim_val <- c(spcov_orig2optim_val, range_log = range_log)
    spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_log = spcov_initial$is_known[["range"]])
  }
  # range <- spcov_initial$initial[["range"]]
  # range_log <- log(range)
  # spcov_orig2optim_val <- c(spcov_orig2optim_val, range_log = range_log)
  # spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_log = spcov_initial$is_known[["range"]])

  # extra p logodds (for now)
  extra <- spcov_initial$initial[["extra"]]
  extra_half <- extra / 2 # because maximum value is 2
  extra_logodds <- logit(extra_half)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, extra_logodds = extra_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, extra_logodds = spcov_initial$is_known[["extra"]])


  # anisotropy parameters
  ## rotate (between 0 and pi radians)
  rotate <- spcov_initial$initial[["rotate"]]
  rotate_prop <- rotate / pi
  rotate_logodds <- logit(rotate_prop)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, rotate_logodds = rotate_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, rotate_logodds = spcov_initial$is_known[["rotate"]])

  ## scale (between 0 and 1)
  scale <- spcov_initial$initial[["scale"]]
  scale_logodds <- logit(scale)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, scale_logodds = scale_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, scale_logodds = spcov_initial$is_known[["scale"]])

  # return covariance parameter vector
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val > 50 & !spcov_orig2optim_is_known, 50, spcov_orig2optim_val)
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val < -50 & !spcov_orig2optim_is_known, -50, spcov_orig2optim_val)

  # return list
  spcov_orig2optim_val <- list(
    value = spcov_orig2optim_val,
    is_known = spcov_orig2optim_is_known,
    n_est = sum(!spcov_orig2optim_is_known)
  )

  # give class to vector
  new_spcov_orig2optim_val <- structure(spcov_orig2optim_val, class = class(spcov_initial))
}

#' @export
spcov_orig2optim.car <- function(spcov_initial, spcov_profiled, data_object, ...) {
  # are variance parameters spcov_profiled
  if (spcov_profiled) {
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
  range <- spcov_initial$initial[["range"]]
  range <- (range - data_object$rho_lb) / (data_object$rho_ub - data_object$rho_lb) # scale to 0,1
  # range <- (range + 1) / 2 # (from -1, 1 to 0, 2 to 0, 1)
  range_logodds <- logit(range)
  spcov_orig2optim_val <- c(spcov_orig2optim_val, range_logodds = range_logodds)
  spcov_orig2optim_is_known <- c(spcov_orig2optim_is_known, range_logodds = spcov_initial$is_known[["range"]])

  # return covariance parameter vector
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val > 50 & !spcov_orig2optim_is_known, 50, spcov_orig2optim_val)
  spcov_orig2optim_val <- ifelse(spcov_orig2optim_val < -50 & !spcov_orig2optim_is_known, -50, spcov_orig2optim_val)

  # return list
  spcov_orig2optim_val <- list(
    value = spcov_orig2optim_val,
    is_known = spcov_orig2optim_is_known,
    n_est = sum(!spcov_orig2optim_is_known)
  )

  # give class to vector
  new_spcov_orig2optim_val <- structure(spcov_orig2optim_val, class = class(spcov_initial))
}
#' @export
spcov_orig2optim.sar <- spcov_orig2optim.car
