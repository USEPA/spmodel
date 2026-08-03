# shared building blocks for the use_gloglik*/use_laploglik* family
# (R/use_gloglik.R, R/use_gloglik_anis.R, R/use_gloglik_known.R,
# R/use_gloglik_known_anis.R, R/use_gloglik_iid.R, R/use_laploglik.R,
# R/use_laploglik_anis.R, R/use_laploglik_known.R,
# R/use_laploglik_known_anis.R) and, where noted, their per-iteration
# objective-function siblings (R/gloglik.R, R/gloglik_anis.R,
# R/laploglik.R, R/laploglik_anis.R). Extracted because these 9+4 files
# duplicate several multi-line blocks verbatim across the {isotropy,
# anisotropy} x {known, estimated} x {Gaussian, Laplace} axes -- see
# R/spcov_transform_helpers.R for the precedent this follows.

#' Build the anisotropy-corrected distance matrix list for known covariance parameters
#'
#' @param spcov_params_val A \code{spcov_params} object with fixed \code{rotate}/\code{scale} values
#' @param data_object The data object
#'
#' @return A list of distance matrices, one per partition, computed after
#'   rotating/rescaling coordinates by the fixed anisotropy parameters
#'
#' @noRd
build_anis_dist_matrix_list <- function(spcov_params_val, data_object) {
  new_coords_list <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
    rotate = spcov_params_val[["rotate"]], scale = spcov_params_val[["scale"]]
  )
  lapply(new_coords_list, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))
}

#' Build a stub optim() output for a fit with no optim() search
#'
#' @param value The -2*log-likelihood (or -2*restricted-log-likelihood) value at
#'   the fixed/known parameters
#'
#' @return A list matching \code{optim()}'s output shape, with every diagnostic
#'   field set to \code{NA} since \code{optim()} was never called
#'
#' @noRd
known_optim_output_stub <- function(value) {
  list(
    method = NA, control = NA, value = value,
    counts = NA, convergence = NA,
    message = NA, hessian = NA
  )
}

#' Trim a converged optim() output down to its reportable fields
#'
#' @param optim_output The raw list returned by \code{optim()}
#' @param optim_dotlist The optim dotlist (supplies \code{method}, \code{control},
#'   and whether \code{hessian} was requested)
#'
#' @return \code{optim_output}, with \code{method}/\code{control} taken from
#'   \code{optim_dotlist} and \code{hessian} dropped to \code{FALSE} unless it was
#'   actually requested
#'
#' @noRd
trim_optim_output <- function(optim_output, optim_dotlist) {
  list(
    method = optim_dotlist$method,
    control = optim_dotlist$control, value = optim_output$value,
    counts = optim_output$counts, convergence = optim_output$convergence,
    message = optim_output$message,
    hessian = if (optim_dotlist$hessian) optim_output$hessian else FALSE
  )
}

#' Unpack an optim-scale parameter vector into spcov/randcov params on their original scale
#'
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param randcov_orig2optim A \code{randcov_orig2optim} object (or \code{NULL})
#' @param par The optim-scale parameter vector (spcov + randcov only -- any
#'   dispersion element must already have been peeled off via
#'   \code{peel_dispersion()})
#' @param spcov_profiled Whether the overall spatial variance is profiled
#' @param randcov_profiled Whether the random effect variances are profiled out
#' @param data_object The data object
#' @param spcov_initial An \code{spcov_initial} object, or \code{NULL}. When
#'   supplied, an estimated \code{ie} that wandered below the numerical floor
#'   \code{spcov_matrix.*()} uses internally is reconciled via
#'   \code{floor_estimated_ie()} so the reported value matches what was
#'   actually used to build Sigma -- only meaningful right after \code{optim()}
#'   converges (the four laploglik-family "estimated" call sites), never
#'   during the per-iteration objective evaluation itself
#'
#' @return A list with \code{spcov_params_val} and \code{randcov_params_val}
#'
#' @noRd
unpack_optim2orig <- function(spcov_orig2optim, randcov_orig2optim, par, spcov_profiled, randcov_profiled,
                              data_object, spcov_initial = NULL) {
  spcov_orig_val <- spcov_optim2orig(spcov_orig2optim, par, spcov_profiled = spcov_profiled, data_object = data_object)
  spcov_params_val <- get_spcov_params(spcov_type = class(spcov_orig2optim), spcov_orig_val = spcov_orig_val)

  if (!is.null(spcov_initial)) {
    spcov_params_val <- floor_estimated_ie(spcov_params_val, spcov_initial$is_known, data_object$diagtol)
  }

  randcov_orig_val <- randcov_optim2orig(randcov_orig2optim, spcov_orig2optim, par,
    randcov_profiled = randcov_profiled,
    spcov_optim2orig = spcov_params_val
  )

  # when the overall spatial variance is profiled out alongside random effect
  # variances, randcov_optim2orig() returns a list bundling both the updated
  # spcov params and the filled-in random effect params -- unpack it here
  if (!is.null(randcov_profiled) && randcov_profiled) {
    spcov_params_val <- randcov_orig_val$spcov_optim2orig
    randcov_orig_val <- randcov_orig_val$fill_orig_val
  }

  randcov_params_val <- randcov_params(randcov_orig_val)

  list(spcov_params_val = spcov_params_val, randcov_params_val = randcov_params_val)
}

#' Peel the dispersion element(s) off the front of a laploglik-family par vector
#'
#' @param dispersion_orig2optim A \code{dispersion_orig2optim} object
#' @param par The current (or converged) optim parameter vector, packed as
#'   \code{[dispersion | spcov | randcov]}
#' @param family The GLM family
#'
#' @return A list with \code{dispersion_params_val} (the dispersion parameters on
#'   their original scale) and \code{par} (the remaining spcov/randcov elements)
#'
#' @noRd
peel_dispersion <- function(dispersion_orig2optim, par, family) {
  dispersion_orig_val <- dispersion_optim2orig(dispersion_orig2optim, par)
  # dispersion_params() inspects substitute(family) to support both unquoted
  # (dispersion_initial(poisson, ...)) and quoted (dispersion_initial("poisson",
  # ...)) call styles -- passing a bare symbol straight through would trip that
  # same unquoted-name branch here and deparse the literal word "family"
  # instead of its value, so it must be wrapped in a call, not passed bare
  dispersion_params_val <- dispersion_params(as.character(family), dispersion_orig_val$fill_orig_val)
  list(dispersion_params_val = dispersion_params_val, par = dispersion_orig_val$new_par)
}

#' Assemble the optim() starting-value vector from spcov/randcov/dispersion pieces
#'
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param randcov_orig2optim A \code{randcov_orig2optim} object, or \code{NULL}
#'   if there are no random effects
#' @param dispersion_orig2optim A \code{dispersion_orig2optim} object, or
#'   \code{NULL} for Gaussian (non-GLM) estimation, which has no dispersion
#'   parameter
#'
#' @return A numeric vector of only the parameters being estimated (i.e.,
#'   excluding any fixed at a known value), packed as
#'   \code{[spcov | randcov | dispersion]} -- dispersion, when present, is
#'   always last, since \code{peel_dispersion()}/\code{dispersion_optim2orig()}
#'   pull it off the end of \code{optim()}'s returned parameter vector
#'
#' @noRd
assemble_optim_par <- function(spcov_orig2optim, randcov_orig2optim = NULL, dispersion_orig2optim = NULL) {
  # subset to !is_known: parameters the user fixed at a known value are held
  # constant rather than passed to optim() as free parameters to search over
  if (is.null(randcov_orig2optim)) {
    par <- spcov_orig2optim$value[!spcov_orig2optim$is_known]
  } else {
    spcov_pars <- spcov_orig2optim$value[!spcov_orig2optim$is_known]
    randcov_pars <- randcov_orig2optim$value[!randcov_orig2optim$is_known]
    par <- c(spcov_pars, randcov_pars)
  }
  if (!is.null(dispersion_orig2optim)) {
    dispersion_pars <- dispersion_orig2optim$value[!dispersion_orig2optim$is_known]
    par <- c(par, dispersion_pars)
  }
  par
}

#' Rescale profiled spatial covariance and random effect variances by the closed-form overall variance
#'
#' @param spcov_params_val A \code{spcov_params} object, on the profiled
#'   (correlation) scale
#' @param randcov_params_val A \code{randcov_params} object, or \code{NULL}
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param dist_matrix_list A list of distance matrices
#' @param spcov_profiled Whether the overall spatial variance is profiled
#' @param randcov_profiled Whether the random effect variances are profiled out
#'
#' @return A list with \code{spcov_params_val} and \code{randcov_params_val},
#'   rescaled onto their final absolute-variance scale. A no-op (returned
#'   unchanged) unless \code{spcov_profiled} and (\code{randcov_profiled} is
#'   \code{NULL} or \code{TRUE}) -- "profiling" means \code{optim()} only
#'   searched over correlation-type parameters (e.g., range) while holding the
#'   overall variance fixed at 1; the actual variance has a closed-form
#'   maximizer given those parameters, computed here and then scaled back into
#'   every variance term
#'
#' @noRd
rescale_profiled_variance <- function(spcov_params_val, randcov_params_val, data_object, estmethod,
                                      dist_matrix_list, spcov_profiled, randcov_profiled) {
  if (spcov_profiled && (is.null(randcov_profiled) ||
    (!is.null(randcov_profiled) && randcov_profiled))) {
    sigma2 <- get_prof_sigma2(
      spcov_params_val, data_object, estmethod,
      dist_matrix_list, randcov_params_val
    )

    spcov_params_val[["de"]] <- sigma2 * spcov_params_val[["de"]]
    spcov_params_val[["ie"]] <- sigma2 * spcov_params_val[["ie"]]

    if (!is.null(randcov_profiled)) {
      randcov_params_val <- sigma2 * randcov_params_val
    }

    if (inherits(spcov_params_val, c("car", "sar"))) {
      spcov_params_val[["extra"]] <- sigma2 * spcov_params_val[["extra"]]
    }
  }

  list(spcov_params_val = spcov_params_val, randcov_params_val = randcov_params_val)
}

#' Evaluate both anisotropy rotation candidates (shared by Layer 1 and Layer 2)
#'
#' Builds both candidate distance matrix lists -- the rotate angle and its
#' pi-complement (\code{abs(pi - rotate)}) -- and evaluates minus-twice the
#' log-likelihood at each, dispatching to the Gaussian or Laplace
#' products/loss functions depending on whether \code{dispersion_params_val}
#' is supplied. This is the shared q1/q2-evaluation primitive underneath both
#' the per-iteration (Layer 1, inside \code{gloglik_anis()}/
#' \code{laploglik_anis()}) and post-hoc (Layer 2, \code{resolve_anis_rotation()})
#' rotation-ambiguity evaluations, extracted so the two layers can't
#' numerically drift apart from each other. See \code{gloglik_anis.R}/
#' \code{laploglik_anis.R} for why this two-candidate evaluation exists at all
#' (a deliberate optimizer-robustness heuristic against getting stuck near the
#' \code{(0, pi)} boundary of the logit-mapped \code{rotate} parameter, not an
#' identifiability correction -- the two candidates generally describe
#' different anisotropy ellipses).
#'
#' @param spcov_params_val A \code{spcov_params} object
#' @param randcov_params_val A \code{randcov_params} object, or \code{NULL}
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param spcov_profiled Whether the overall spatial variance is profiled
#' @param randcov_profiled Whether the random effect variances are profiled out
#' @param dispersion_params_val A \code{dispersion_params} object, or
#'   \code{NULL} for Gaussian (non-GLM) estimation -- when supplied, dispatches
#'   to the Laplace-approximation products/loss functions instead of the
#'   Gaussian ones
#'
#' @return A list with \code{dist_matrix_list_q1}, \code{dist_matrix_list_q2},
#'   \code{minustwo_q1}, and \code{minustwo_q2}
#'
#' @noRd
anis_rotation_candidates <- function(spcov_params_val, randcov_params_val, data_object, estmethod,
                                     spcov_profiled, randcov_profiled, dispersion_params_val = NULL) {
  new_coords_list_q1 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
    rotate = spcov_params_val[["rotate"]], scale = spcov_params_val[["scale"]]
  )
  dist_matrix_list_q1 <- lapply(new_coords_list_q1, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

  new_coords_list_q2 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
    rotate = abs(pi - spcov_params_val[["rotate"]]), scale = spcov_params_val[["scale"]]
  )
  dist_matrix_list_q2 <- lapply(new_coords_list_q2, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

  if (is.null(dispersion_params_val)) {
    prods_q1 <- gloglik_products(spcov_params_val, data_object, estmethod, dist_matrix_list_q1, randcov_params_val)
    prods_q2 <- gloglik_products(spcov_params_val, data_object, estmethod, dist_matrix_list_q2, randcov_params_val)
    minustwo_q1 <- get_minustwologlik(prods_q1, estmethod, data_object$n, data_object$p, spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled)
    minustwo_q2 <- get_minustwologlik(prods_q2, estmethod, data_object$n, data_object$p, spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled)
  } else {
    prods_q1 <- laploglik_products(spcov_params_val, dispersion_params_val, data_object, estmethod, dist_matrix_list_q1, randcov_params_val)
    prods_q2 <- laploglik_products(spcov_params_val, dispersion_params_val, data_object, estmethod, dist_matrix_list_q2, randcov_params_val)
    minustwo_q1 <- get_minustwolaploglik(prods_q1, estmethod, data_object$n, data_object$p, spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled)
    minustwo_q2 <- get_minustwolaploglik(prods_q2, estmethod, data_object$n, data_object$p, spcov_profiled = spcov_profiled, randcov_profiled = randcov_profiled)
  }

  list(
    dist_matrix_list_q1 = dist_matrix_list_q1, dist_matrix_list_q2 = dist_matrix_list_q2,
    minustwo_q1 = minustwo_q1, minustwo_q2 = minustwo_q2
  )
}

#' Resolve the post-hoc (Layer 2) anisotropy rotation ambiguity after optim() converges
#'
#' gloglik_anis()/laploglik_anis() (the per-iteration objective optim() calls)
#' already evaluate both the rotate angle and its pi-complement at every step
#' and report whichever fits better -- see \code{anis_rotation_candidates()}
#' for the shared evaluation and why it exists. Now that optim() has
#' converged, this redoes that same two-candidate evaluation once more, but
#' this time keeps track of *which* candidate wins so the reported
#' spcov_params_val and dist_matrix_list can be permanently fixed to match.
#' This is Layer 2 only -- the per-iteration Layer 1 evaluation inside
#' gloglik_anis()/laploglik_anis() is a separate call site with a different
#' return contract (scalar min only vs quadrant-index-plus-params) and must
#' keep evaluating both candidates at every iteration.
#'
#' @param spcov_params_val A \code{spcov_params} object
#' @param randcov_params_val A \code{randcov_params} object, or \code{NULL}
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param spcov_profiled Whether the overall spatial variance is profiled
#' @param randcov_profiled Whether the random effect variances are profiled out
#' @param dispersion_params_val A \code{dispersion_params} object, or
#'   \code{NULL} for Gaussian (non-GLM) estimation -- when supplied, dispatches
#'   to the Laplace-approximation products/loss functions instead of the
#'   Gaussian ones
#'
#' @return A list with the (possibly rotate-flipped) \code{spcov_params_val}
#'   and the winning candidate's \code{dist_matrix_list}
#'
#' @noRd
resolve_anis_rotation <- function(spcov_params_val, randcov_params_val, data_object, estmethod,
                                  spcov_profiled, randcov_profiled, dispersion_params_val = NULL) {
  candidates <- anis_rotation_candidates(
    spcov_params_val, randcov_params_val, data_object, estmethod,
    spcov_profiled, randcov_profiled, dispersion_params_val
  )

  # lower -2loglik (i.e., higher likelihood) wins between the two candidates
  rotate_min <- which.min(c(candidates$minustwo_q1, candidates$minustwo_q2))

  if (rotate_min == 1) {
    dist_matrix_list <- candidates$dist_matrix_list_q1
  } else if (rotate_min == 2) {
    spcov_params_val[["rotate"]] <- abs(pi - spcov_params_val[["rotate"]])
    dist_matrix_list <- candidates$dist_matrix_list_q2
  }

  list(spcov_params_val = spcov_params_val, dist_matrix_list = dist_matrix_list)
}
