#' Compute small inverse via Helmert-Wolf blocking
#'
#' @param SigInv An inverse covariance matrix
#' @param Sigldet A log determinant
#' @param observed_index Index of observed values
#'
#' @return A small inverse via Helmert-Wolf blocking
#'
#' @noRd
# SigInv/Sigldet are computed for the full (observed + missing) precision
# matrix; when some locations are unobserved, this collapses them down to the
# inverse and log determinant of the marginal covariance for only the
# observed locations, via a Schur-complement (Helmert-Wolf) block reduction
# instead of inverting the full matrix and re-extracting a submatrix
hwInv <- function(SigInv, Sigldet, observed_index = NULL) {
  # only need to do the block reduction if there actually are missing locations
  if (NROW(SigInv) > length(observed_index)) {
    missing_index <- which(!(seq_len(NROW(SigInv)) %in% observed_index))
    # partition the joint precision matrix into observed/missing blocks
    SigInv_oo <- SigInv[observed_index, observed_index, drop = FALSE]
    SigInv_om <- SigInv[observed_index, missing_index, drop = FALSE]
    SigInv_mm <- SigInv[missing_index, missing_index, drop = FALSE]
    SigInv_mm_upchol <- chol(forceSymmetric(SigInv_mm))
    # log determinant of the full matrix = log det(SigInv_mm) + log det(Schur complement),
    # so subtracting out the missing-block contribution leaves the observed-only log determinant
    Sigldet <- Sigldet + 2 * sum(log(diag(SigInv_mm_upchol)))
    # Schur complement: SigInv_oo - SigInv_om %*% SigInv_mm^(-1) %*% t(SigInv_om)
    # is the precision matrix of the observed locations alone
    SigInv <- SigInv_oo - SigInv_om %*% tcrossprod(chol2inv(SigInv_mm_upchol), SigInv_om)
  }
  list(SigInv = SigInv, Sigldet = Sigldet)
}
