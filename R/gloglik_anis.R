#' Find (minus twice negative) Gaussian log-likelihood while optimizing (with anisotropy)
#'
#' @param par Parameters to optimize over
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param estmethod The estimation method
#' @param X Model matrix
#' @param y Response vector
#' @param n Sample size
#' @param p Number of fixed effects
#' @param xcoord_val A vector of x-coordinates
#' @param ycoord_val A vector of y-coordinates
#' @param spcov_profiled Whether the overall spatial variance is profiled
#' @param randcov_orig2optim A \code{randcov_orig2optim} object
#' @param randcov_Zs Random effect design matrices
#' @param observed_index Index of observed values
#' @param partition_matrix Partition matrix
#'
#' @return (Minus twice negative) Gaussian log-likelihood (with anisotropy)
#'
#' @noRd
gloglik_anis <- function(par, spcov_orig2optim, data_object, estmethod,
                         spcov_profiled, randcov_orig2optim = NULL,
                         randcov_profiled = NULL) {
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
  candidates <- anis_rotation_candidates(spcov_params_val, randcov_params_val, data_object, estmethod, spcov_profiled, randcov_profiled)

  # keep whichever rotation quadrant gives the better (smaller) -2loglik
  minustwologlik <- min(c(candidates$minustwo_q1, candidates$minustwo_q2))
}
