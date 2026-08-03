#' Evaluate minus twice the Laplace-approximated log-likelihood under anisotropy
#'
#' @param par The current optimization parameter vector
#' @param spcov_orig2optim A spatial covariance parameter list on the optimization scale
#' @param dispersion_orig2optim A dispersion parameter list on the optimization scale
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param spcov_profiled Whether the spatial covariance parameters are profiled out
#' @param randcov_orig2optim A random effect variance list on the optimization
#'   scale (or \code{NULL} if there are no random effects)
#' @param randcov_profiled Whether the random effect variances are profiled out
#'
#' @return The smaller of the two minus-twice log-likelihood values obtained by
#'   evaluating the rotation angle and its `pi`-complement -- an
#'   optimizer-robustness heuristic (see the inline comment below), not an
#'   identifiability correction; the two candidates generally give different
#'   values
#'
#' @noRd
laploglik_anis <- function(par, spcov_orig2optim, dispersion_orig2optim, data_object, estmethod,
                           spcov_profiled, randcov_orig2optim = NULL,
                           randcov_profiled = NULL) {
  # optim() passes a single flat parameter vector; the dispersion element(s)
  # are pulled off the front first and converted back to their original
  # scale, leaving the remaining (spatial/random effect) parameters in par
  # dispersion first then remove
  dispersion_result <- peel_dispersion(dispersion_orig2optim, par, data_object$family)
  dispersion_params_val <- dispersion_result$dispersion_params_val
  par <- dispersion_result$par

  unpacked <- unpack_optim2orig(spcov_orig2optim, randcov_orig2optim, par, spcov_profiled, randcov_profiled, data_object)
  spcov_params_val <- unpacked$spcov_params_val
  randcov_params_val <- unpacked$randcov_params_val

  # evaluate the objective at both the estimated rotate angle and its
  # pi-complement (abs(pi - rotate)) and keep whichever fits better. These two
  # candidates generally describe different anisotropy ellipses -- this is a
  # deliberate optimizer-robustness heuristic, not an identifiability
  # correction. rotate is mapped from an unconstrained optim() scale via a
  # logit-type transform that flattens sharply near the (0, pi) boundary, so a
  # single-candidate search can get stuck making vanishingly slow progress
  # near either edge. Reporting min(g(rotate), g(pi-rotate)) lets optim()'s
  # own rotate value converge anywhere in (0, pi): if the better-fitting
  # ellipse is actually near the opposite edge, its mirror candidate is
  # evaluated and rewarded throughout the search, so optim() never has to
  # cross the flat region to find it. See anis_rotation_candidates() in
  # R/use_loglik_helpers.R for the shared two-candidate evaluation (also used
  # by the post-hoc Layer 2 resolution, resolve_anis_rotation()).
  candidates <- anis_rotation_candidates(
    spcov_params_val, randcov_params_val, data_object, estmethod,
    spcov_profiled, randcov_profiled, dispersion_params_val
  )

  # keep whichever rotation quadrant gives the better (smaller) -2loglik
  minustwolaploglik <- min(c(candidates$minustwo_q1, candidates$minustwo_q2))
}
