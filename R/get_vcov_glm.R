#' Assemble the fixed effect variance-covariance list for a GLM-type model
#'
#' @param cov_betahat_corrected The covariance matrix of betahat, corrected
#'   for estimation of the latent random effects (see \code{get_wts_varw()})
#' @param cov_betahat_uncorrected The uncorrected covariance matrix of betahat
#'
#' @return A list with element \code{fixed}, itself a list with elements
#'   \code{corrected} and \code{uncorrected}
#'
#' @noRd
get_vcov_glm <- function(cov_betahat_corrected, cov_betahat_uncorrected) {
  # negative variances (diagonal entries) signal numerical instability in the
  # fit, often from a near-singular covariance matrix with tiny independent error
  if (any(diag(cov_betahat_corrected) < 0)) {
    warning("Model fit potentially unstable. Consider fixing ie (via spcov_initial) at some non-zero value greater than 1e-4 and refitting the model.", call. = FALSE)
  }
  vcov_fixed <- list(corrected = cov_betahat_corrected, uncorrected = cov_betahat_uncorrected)
  list(fixed = vcov_fixed)
}
