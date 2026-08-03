#' Compute the inverse covariance sum with random effects
#'
#' @param Ainv The inverse of a matrix A
#' @param Aldet The log determinant of a matrix A
#' @param randcov_params A \code{randcov_params} object
#' @param randcov_Zs Random effect design matrices
#'
#' @return The inverse of A + B, where B is the random effects covariance
#' (only used for autoregressive models, which are parameterized in terms
#' of the precision matrix).
#'
#' @noRd
smwInv_rand <- function(Ainv, Aldet, randcov_params, randcov_Zs) {
  if (!is.null(randcov_params)) {
    # each random effect adds a low-rank term Z %*% (var * I) %*% t(Z) to the
    # covariance, so its inverse and log-determinant are folded in one term at
    # a time via the Sherman-Morrison-Woodbury (SMW) identity below, rather than
    # forming and inverting the full (dense) covariance matrix directly
    for (i in seq_along(randcov_params)) {
      randcov_var <- names(randcov_params)[[i]]
      Z <- randcov_Zs[[randcov_var]][["Z"]]
      Ainv_Z <- Ainv %*% Z
      smw_mid <- crossprod(Z, Ainv_Z)
      # smw_mid = (1 / var) * I + t(Z) %*% Ainv %*% Z -- the "middle" matrix
      # inverted by the SMW identity, whose dimension is the (small) number of
      # random effect levels rather than the (large) number of observations
      diag(smw_mid) <- 1 / randcov_params[[randcov_var]] + diag(smw_mid)
      smw_mid_upchol <- chol(forceSymmetric(smw_mid))
      Inv_smw_mid <- chol2inv(smw_mid_upchol)
      # log-det of smw_mid from its Cholesky factor (2 * sum log diag), reused
      # below via the matrix determinant lemma
      ldet_smw_mid <- 2 * sum(log(diag(smw_mid_upchol)))
      # SMW inverse update: (A + Z var Z')^-1 = Ainv - Ainv Z (smw_mid)^-1 Z' Ainv
      Ainv <- Ainv - tcrossprod(Ainv_Z %*% Inv_smw_mid, Ainv_Z)
      # matrix determinant lemma: log|A + Z var Z'| = log|A| + log|var| * ncol(Z) + log|smw_mid|
      Aldet <- Aldet + NCOL(Z) * log(randcov_params[[randcov_var]]) + ldet_smw_mid
    }
  }
  list(SigInv = Ainv, Sigldet = Aldet)
}
