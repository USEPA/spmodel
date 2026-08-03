# shared building blocks for cov_estimate_gloglik.R / cov_estimate_laploglik.R.
# Each of the 4 dispatcher functions there (cov_estimate_gloglik_splm/
# _spautor, cov_estimate_laploglik_spglm/_spgautor) used to write its own
# "which use_gloglik*/use_laploglik* variant to call" decision tree twice --
# once for models with no random effects, once for models with random
# effects -- differing only in whether randcov_initial_val is NULL. Since
# all()/any() in R silently drop NULL arguments, a single decision tree
# parameterized by randcov_initial_val (NULL-safe) reproduces both branches
# exactly; splm/spautor and glm/gautor keep separate helpers below because
# their case logic (branch order, extra known-ness checks) genuinely differs,
# not just their random-effects handling.
#
# IMPORTANT: randcov_profiled must be NULL (not FALSE) whenever
# randcov_initial_val is NULL, not just "roughly false" -- inside use_gloglik()/
# use_gloglik_anis(), the profiled-variance rescale is gated by
# `spcov_profiled && (is.null(randcov_profiled) || (!is.null(randcov_profiled)
# && randcov_profiled))`, which evaluates differently for NULL vs FALSE when
# spcov_profiled is TRUE. Verified against every existing call site: whenever
# randcov_profiled was previously passed explicitly (not omitted), its value
# always equalled spcov_profiled -- hence the `if (is.null(randcov_initial_val))
# NULL else spcov_profiled` rule used throughout below.

#' Choose and run the Gaussian log-likelihood estimator for a geostatistical model
#'
#' @param spcov_initial_val An \code{spcov_initial} object with initial values filled in
#' @param randcov_initial_val A \code{randcov_initial} object, or \code{NULL} if no random effects
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param dist_matrix_list A distance matrix list (or \code{NULL} under anisotropy)
#' @param optim_dotlist The optim dotlist
#'
#' @return The Gaussian log-likelihood estimates (mirrors \code{cov_estimate_gloglik_splm()}'s
#'   prior inline branching)
#'
#' @noRd
run_gloglik_dispatch_splm <- function(spcov_initial_val, randcov_initial_val, data_object, estmethod,
                                      dist_matrix_list, optim_dotlist) {
  # the iid shortcut is only ever checked when there's no random effect --
  # matches that the prior randcov-present branch had no iid case at all
  allow_iid <- is.null(randcov_initial_val)

  de_known <- spcov_initial_val$is_known[["de"]]
  de_known_zero <- de_known && (spcov_initial_val$initial[["de"]] == 0)
  ie_known <- spcov_initial_val$is_known[["ie"]]

  if (all(spcov_initial_val$is_known, randcov_initial_val$is_known)) {
    if (data_object$anisotropy) {
      use_gloglik_known_anis(spcov_initial_val, data_object, estmethod, randcov_initial_val)
    } else {
      use_gloglik_known(spcov_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial_val)
    }
  } else if (allow_iid && de_known_zero) {
    use_gloglik_iid(spcov_initial_val, estmethod, data_object, dist_matrix_list)
  } else if (any(de_known && !de_known_zero, ie_known, randcov_initial_val$is_known)) {
    run_gloglik_optim_splm(spcov_initial_val, randcov_initial_val, data_object, estmethod, dist_matrix_list, optim_dotlist, spcov_profiled = FALSE)
  } else {
    run_gloglik_optim_splm(spcov_initial_val, randcov_initial_val, data_object, estmethod, dist_matrix_list, optim_dotlist, spcov_profiled = TRUE)
  }
}

#' Run \code{use_gloglik()}/\code{use_gloglik_anis()} with the correct randcov_profiled
#'
#' @param spcov_profiled Is spatial profiling used for this call?
#' @inheritParams run_gloglik_dispatch_splm
#'
#' @noRd
run_gloglik_optim_splm <- function(spcov_initial_val, randcov_initial_val, data_object, estmethod,
                                   dist_matrix_list, optim_dotlist, spcov_profiled) {
  randcov_profiled_val <- if (is.null(randcov_initial_val)) NULL else spcov_profiled
  if (data_object$anisotropy) {
    use_gloglik_anis(spcov_initial_val, data_object, estmethod,
      spcov_profiled = spcov_profiled,
      randcov_initial = randcov_initial_val, randcov_profiled = randcov_profiled_val,
      optim_dotlist = optim_dotlist
    )
  } else {
    use_gloglik(spcov_initial_val, data_object, estmethod, dist_matrix_list,
      spcov_profiled = spcov_profiled,
      randcov_initial = randcov_initial_val, randcov_profiled = randcov_profiled_val,
      optim_dotlist = optim_dotlist
    )
  }
}

#' Choose and run the Gaussian log-likelihood estimator for an autoregressive model
#'
#' @inheritParams run_gloglik_dispatch_splm
#'
#' @return The Gaussian log-likelihood estimates (mirrors \code{cov_estimate_gloglik_spautor()}'s
#'   prior inline branching -- note the branch order differs from the splm version:
#'   car/sar's third variance component (\code{extra}) means the profiled case is
#'   checked before the non-profiled case, the opposite order from splm)
#'
#' @noRd
run_gloglik_dispatch_spautor <- function(spcov_initial_val, randcov_initial_val, data_object, estmethod,
                                         dist_matrix_list, optim_dotlist) {
  allow_iid <- is.null(randcov_initial_val)

  de_known <- spcov_initial_val$is_known[["de"]]
  de_known_zero <- de_known && (spcov_initial_val$initial[["de"]] == 0)
  ie_known <- spcov_initial_val$is_known[["ie"]]
  extra_known <- spcov_initial_val$is_known[["extra"]]
  extra_known_zero <- extra_known && (spcov_initial_val$initial[["extra"]] == 0)

  if (all(spcov_initial_val$is_known, randcov_initial_val$is_known)) {
    use_gloglik_known(spcov_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial_val)
  } else if (allow_iid && de_known_zero && extra_known_zero) {
    use_gloglik_iid(spcov_initial_val, estmethod, data_object, dist_matrix_list)
  } else if (extra_known_zero && !any(de_known, ie_known, randcov_initial_val$is_known)) {
    run_gloglik_optim_spautor(spcov_initial_val, randcov_initial_val, data_object, estmethod, dist_matrix_list, optim_dotlist, spcov_profiled = TRUE)
  } else {
    run_gloglik_optim_spautor(spcov_initial_val, randcov_initial_val, data_object, estmethod, dist_matrix_list, optim_dotlist, spcov_profiled = FALSE)
  }
}

#' Run \code{use_gloglik()} (spautor never has anisotropy) with the correct randcov_profiled
#'
#' @param spcov_profiled Is spatial profiling used for this call?
#' @inheritParams run_gloglik_dispatch_splm
#'
#' @noRd
run_gloglik_optim_spautor <- function(spcov_initial_val, randcov_initial_val, data_object, estmethod,
                                      dist_matrix_list, optim_dotlist, spcov_profiled) {
  randcov_profiled_val <- if (is.null(randcov_initial_val)) NULL else spcov_profiled
  use_gloglik(spcov_initial_val, data_object, estmethod, dist_matrix_list,
    spcov_profiled = spcov_profiled,
    randcov_initial = randcov_initial_val, randcov_profiled = randcov_profiled_val,
    optim_dotlist = optim_dotlist
  )
}

#' Choose and run the Laplace log-likelihood estimator for a GLM model
#'
#' @param spcov_initial_val An \code{spcov_initial} object with initial values filled in
#' @param dispersion_initial_val A \code{dispersion_initial} object with initial values filled in
#' @param randcov_initial_val A \code{randcov_initial} object, or \code{NULL} if no random effects
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param dist_matrix_list A distance matrix list (or \code{NULL} under anisotropy)
#' @param optim_dotlist The optim dotlist
#'
#' @return The Laplace log-likelihood estimates (mirrors \code{cov_estimate_laploglik_spglm()}'s
#'   prior inline branching -- unlike the Gaussian family, there is no iid shortcut and no
#'   profiled variant, since \code{laploglik()} never supports profiling)
#'
#' @noRd
run_laploglik_dispatch_spglm <- function(spcov_initial_val, dispersion_initial_val, randcov_initial_val,
                                         data_object, estmethod, dist_matrix_list, optim_dotlist) {
  de_known <- spcov_initial_val$is_known[["de"]]
  ie_known <- spcov_initial_val$is_known[["ie"]]
  dispersion_known <- dispersion_initial_val$is_known[["dispersion"]]

  if (all(de_known, ie_known, dispersion_known, randcov_initial_val$is_known)) {
    if (data_object$anisotropy) {
      use_laploglik_known_anis(spcov_initial_val, dispersion_initial_val, data_object, estmethod, randcov_initial = randcov_initial_val)
    } else {
      use_laploglik_known(spcov_initial_val, dispersion_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial_val)
    }
  } else {
    # spcov_profiled is always FALSE for the Laplace family, so the NULL-vs-
    # FALSE rescale-gate subtlety documented at the top of this file never
    # actually triggers here -- randcov_profiled still follows the same rule
    # for consistency with the Gaussian family
    randcov_profiled_val <- if (is.null(randcov_initial_val)) NULL else FALSE
    if (data_object$anisotropy) {
      use_laploglik_anis(spcov_initial_val, dispersion_initial_val, data_object, estmethod,
        spcov_profiled = FALSE,
        randcov_initial = randcov_initial_val, randcov_profiled = randcov_profiled_val,
        optim_dotlist = optim_dotlist
      )
    } else {
      use_laploglik(spcov_initial_val, dispersion_initial_val, data_object, estmethod, dist_matrix_list,
        spcov_profiled = FALSE,
        randcov_initial = randcov_initial_val, randcov_profiled = randcov_profiled_val,
        optim_dotlist = optim_dotlist
      )
    }
  }
}

#' Choose and run the composite-likelihood estimator (Curriero and Lele, 1999)
#'
#' @param spcov_initial_val An \code{spcov_initial} object with initial values filled in
#' @param data_object The data object
#' @param dist_matrix_list A distance matrix list (or \code{NULL} under anisotropy)
#' @param optim_dotlist The optim dotlist
#'
#' @return The composite log-likelihood estimates (mirrors \code{cov_estimate_cl()}'s
#'   prior inline branching)
#'
#' @noRd
run_cl_dispatch <- function(spcov_initial_val, data_object, dist_matrix_list, optim_dotlist) {
  if (all(spcov_initial_val$is_known)) {
    # nothing to estimate -- just evaluate/return the fixed parameters
    use_glogclik_known(spcov_initial_val, data_object, dist_matrix_list, data_object$partition_list)
  } else {
    # optimize the composite log-likelihood over the unknown parameters
    use_glogclik(spcov_initial_val, data_object, dist_matrix_list, data_object$partition_list, optim_dotlist)
  }
}

#' Choose and run the semivariogram weighted-least-squares estimator
#'
#' @param spcov_initial_val An \code{spcov_initial} object with initial values filled in
#' @param data_object The data object
#' @param dist_matrix_list A distance matrix list (or \code{NULL} under anisotropy)
#' @param esv The empirical semivariogram
#' @param weights sv-wls weights
#' @param optim_dotlist The optim dotlist
#'
#' @return Covariance parameter estimates (mirrors \code{cov_estimate_sv()}'s
#'   prior inline branching)
#'
#' @noRd
run_sv_dispatch <- function(spcov_initial_val, data_object, dist_matrix_list, esv, weights, optim_dotlist) {
  if (all(spcov_initial_val$is_known)) {
    # nothing to estimate -- just evaluate the loss at the fixed parameters
    use_svloss_known(spcov_initial_val, dist_matrix_list, esv, weights)
  } else {
    # minimize the weighted least-squares loss over the unknown parameters
    use_svloss(spcov_initial_val, dist_matrix_list, esv, weights, optim_dotlist, data_object = data_object)
  }
}

#' Choose and run the Laplace log-likelihood estimator for an areal (autoregressive) GLM model
#'
#' @inheritParams run_laploglik_dispatch_spglm
#'
#' @return The Laplace log-likelihood estimates (mirrors \code{cov_estimate_laploglik_spgautor()}'s
#'   prior inline branching -- never branches on anisotropy, since areal models have none, and
#'   includes \code{extra} in the known-ness checks for car/sar's third variance component)
#'
#' @noRd
run_laploglik_dispatch_spgautor <- function(spcov_initial_val, dispersion_initial_val, randcov_initial_val,
                                            data_object, estmethod, dist_matrix_list, optim_dotlist) {
  de_known <- spcov_initial_val$is_known[["de"]]
  ie_known <- spcov_initial_val$is_known[["ie"]]
  extra_known <- spcov_initial_val$is_known[["extra"]]
  dispersion_known <- dispersion_initial_val$is_known[["dispersion"]]

  if (all(de_known, ie_known, extra_known, dispersion_known, randcov_initial_val$is_known)) {
    use_laploglik_known(spcov_initial_val, dispersion_initial_val, data_object, estmethod, dist_matrix_list, randcov_initial_val)
  } else {
    randcov_profiled_val <- if (is.null(randcov_initial_val)) NULL else FALSE
    use_laploglik(spcov_initial_val, dispersion_initial_val, data_object, estmethod, dist_matrix_list,
      spcov_profiled = FALSE,
      randcov_initial = randcov_initial_val, randcov_profiled = randcov_profiled_val,
      optim_dotlist = optim_dotlist
    )
  }
}
