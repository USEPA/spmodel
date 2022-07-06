#' Compute small inverse via Helmert-Wolf blocking
#'
#' @param SigInv An inverse covariance matrix
#' @param Sigldet A log determinant
#' @param observed_index Index of observed values
#'
#' @return A small inverse via Helmert-Wolf blocking
#'
#' @noRd
hwInv <- function(SigInv, Sigldet, observed_index = NULL) {
  if (NROW(SigInv) > length(observed_index)) {
    missing_index <- which(!(seq_len(NROW(SigInv)) %in% observed_index))
    SigInv_oo <- SigInv[observed_index, observed_index, drop = FALSE]
    SigInv_om <- SigInv[observed_index, missing_index, drop = FALSE]
    SigInv_mm <- SigInv[missing_index, missing_index, drop = FALSE]
    SigInv_mm_upchol <- chol(forceSymmetric(SigInv_mm))
    Sigldet <- Sigldet + 2 * sum(log(diag(SigInv_mm_upchol)))
    SigInv <- SigInv_oo - SigInv_om %*% tcrossprod(chol2inv(SigInv_mm_upchol), SigInv_om)
  }
  list(SigInv = SigInv, Sigldet = Sigldet)
}
