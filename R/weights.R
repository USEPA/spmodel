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

use_cressie_dr_weights <- function(np, sv_val) {
  np / sv_val
}

use_cressie_nopairs_weights <- function(np, sv_val) {
  1 / sv_val^2
}

use_cressie_dr_nopairs_weights <- function(np, sv_val) {
  1 / sv_val
}

use_pairs_weights <- function(np) {
  np
}

use_pairs_invd_weights <- function(np, dist_vector) {
  np / dist_vector
}

use_pairs_invsd_weights <- function(np, dist_vector) {
  np / dist_vector^2
}
