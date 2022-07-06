#' Get semivariogram loss
#'
#' @param spcov_params A \code{spcov_params} object
#' @param gamma Semivariogram value
#' @param weights wls weights
#' @param dist_vector A distance vector
#' @param np The number of pairs
#'
#' @return The semivariogram loss
#'
#' @noRd
get_svloss <- function(spcov_params, esv, weights) {


  # define esv values
  gamma <- esv[["gamma"]]
  dist_vector <- esv[["dist"]]
  np <- esv[["np"]]

  # rest of function
  sigma2_val <- spcov_params[["de"]] + spcov_params[["ie"]]
  spcov_vec_val <- spcov_vector(spcov_params, dist_vector)
  sv_val <- sigma2_val - spcov_vec_val
  weights_val <- switch(weights,
    "cressie" = use_cressie_weights(np = np, sv_val = sv_val),
    "cressie-dr" = use_cressie_dr_weights(np = np, sv_val = sv_val),
    "cressie-nopairs" = use_cressie_nopairs_weights(np = np, sv_val = sv_val),
    "cressie-dr-nopairs" = use_cressie_dr_nopairs_weights(np = np, sv_val = sv_val),
    "pairs" = use_pairs_weights(np = np),
    "pairs-invd" = use_pairs_invd_weights(np = np, dist_vector = dist_vector),
    "pairs-invrd" = use_pairs_invsd_weights(np = np, dist_vector = dist_vector),
    "ols" = 1,
  )
  wls_val <- sum(weights_val * (gamma - sv_val)^2)
  wls_val
}
