#' Switch functions for sv-wls weights
#'
#' @param np Empirical semivariogram pairs (corresponding to a distance vector)
#' @param sv_val Empirical semivariogram value (corresponding to a distance vector)
#'
#' @return sv-wls weights
#'
#' @noRd
use_cressie_weights <- function(np, sv_val) {
  np / sv_val^2
}

#' Cressie's "directly reweighted" sv-wls weights
#'
#' @param np Empirical semivariogram pairs (corresponding to a distance vector)
#' @param sv_val Empirical semivariogram value (corresponding to a distance vector)
#'
#' @return sv-wls weights
#'
#' @noRd
use_cressie_dr_weights <- function(np, sv_val) {
  np / sv_val
}

#' Cressie's weights, ignoring pair counts
#'
#' @param np Empirical semivariogram pairs (unused; kept for a consistent signature)
#' @param sv_val Empirical semivariogram value (corresponding to a distance vector)
#'
#' @return sv-wls weights
#'
#' @noRd
use_cressie_nopairs_weights <- function(np, sv_val) {
  1 / sv_val^2
}

#' Cressie's "directly reweighted" weights, ignoring pair counts
#'
#' @param np Empirical semivariogram pairs (unused; kept for a consistent signature)
#' @param sv_val Empirical semivariogram value (corresponding to a distance vector)
#'
#' @return sv-wls weights
#'
#' @noRd
use_cressie_dr_nopairs_weights <- function(np, sv_val) {
  1 / sv_val
}

#' Pair-count sv-wls weights
#'
#' @param np Empirical semivariogram pairs (corresponding to a distance vector)
#'
#' @return sv-wls weights
#'
#' @noRd
use_pairs_weights <- function(np) {
  np
}

#' Pair-count weights, scaled by inverse distance
#'
#' @param np Empirical semivariogram pairs (corresponding to a distance vector)
#' @param dist_vector A distance vector
#'
#' @return sv-wls weights
#'
#' @noRd
use_pairs_invd_weights <- function(np, dist_vector) {
  np / dist_vector
}

#' Pair-count weights, scaled by inverse squared distance
#'
#' @param np Empirical semivariogram pairs (corresponding to a distance vector)
#' @param dist_vector A distance vector
#'
#' @return sv-wls weights
#'
#' @noRd
use_pairs_invsd_weights <- function(np, dist_vector) {
  np / dist_vector^2
}
