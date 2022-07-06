#' Compute the inverse covariance sum with random effects
#'
#' @param Ainv The inverse of a matrix A
#' @param Aldet The log determinant of a matrix A
#' @param randcov_params A \code{randcov_params} object
#' @param randcov_Zs Random effect design matrices
#'
#' @return The inverse of A + B, where B is the random effects covariance
#'
#' @noRd
smwInv_rand <- function(Ainv, Aldet, randcov_params, randcov_Zs) {
  if (!is.null(randcov_params)) {
    for (i in seq_along(randcov_params)) {
      randcov_var <- names(randcov_params)[[i]]
      Z <- randcov_Zs[[randcov_var]][["Z"]]
      Ainv_Z <- Ainv %*% Z
      smw_mid <- crossprod(Z, Ainv_Z)
      diag(smw_mid) <- 1 / randcov_params[[randcov_var]] + diag(smw_mid)
      smw_mid_upchol <- chol(forceSymmetric(smw_mid))
      Inv_smw_mid <- chol2inv(smw_mid_upchol)
      ldet_smw_mid <- 2 * sum(log(diag(smw_mid_upchol)))
      Ainv <- Ainv - tcrossprod(Ainv_Z %*% Inv_smw_mid, Ainv_Z)
      Aldet <- Aldet + NCOL(Z) * log(randcov_params[[randcov_var]]) + ldet_smw_mid
    }
  }
  list(SigInv = Ainv, Sigldet = Aldet)
}
