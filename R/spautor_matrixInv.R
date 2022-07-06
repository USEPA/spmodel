# compute matrix inverse for autoregressive models
spautor_cov_matrixInv <- function(spcov_params_val, data_object,
                                  dist_matrix_list, randcov_params_val, ldet = TRUE) {
  if (spcov_params_val[["de"]] < 0.001) {
    spcov_params_val[["de"]] <- 0.001
  }

  # MAKE IT CLEAR DIST_MATRIX_LIST IS NOT A LIST
  dist_matrix <- dist_matrix_list

  if (is.null(data_object$partition_matrix)) {
    # compute the dependent inverse portion
    SigInv_de <- spcov_matrixInv_de(spcov_params_val, dist_matrix, data_object$M)
    # using Matrix if loaded
    # find its cholesky
    # chol_SigInv_de <- chol(forceSymmetric(SigInv_de))
    # ldet_Sig_de <- - 2 * sum(log(diag(chol_SigInv_de)))
    ldet_Sig_de <- -as.numeric(Matrix::determinant(SigInv_de)$modulus)
    # if there is no independent random error variance, stop with the inverse
    # and log determinant -- else find inverse and log determinant of
    # (de + ie)
    if (spcov_params_val[["ie"]] == 0) {
      # if the de is zero this is the inverse
      SigInv <- SigInv_de
      # if the de is zero this is the ldet
      ldet_Sig <- ldet_Sig_de
    } else {
      # finding the cholesky and inverse of middle smw (A^-1 + B^-1)^-1
      ## FOR RAND EFFECTS NEED TO SMW STEP THROUGH Z %*% I %*% Zt
      smw_mid <- SigInv_de
      diag(smw_mid) <- 1 / spcov_params_val[["ie"]] + diag(smw_mid)
      smw_mid_upchol <- chol(forceSymmetric(smw_mid))
      Inv_smw_mid <- chol2inv(smw_mid_upchol)
      ldet_smw_mid <- 2 * sum(log(diag(smw_mid_upchol)))
      # finding the inverse (A + B)^-1 = A^-1 - A^-1 (A^-1 + B^-1)^-1 A^-1
      SigInv <- SigInv_de - SigInv_de %*% Inv_smw_mid %*% SigInv_de
      # finding sample size
      # n <- length(y)
      # finding ldet (A + B)^-1 = A^-1 - A^-1 (A^-1 + B^-1)^-1 A^-1
      ldet_Sig <- ldet_Sig_de + data_object$n * log(spcov_params_val[["ie"]]) + ldet_smw_mid
    }

    # browser()
    # do random effects
    smwInv_rand_val <- smwInv_rand(SigInv, ldet_Sig, randcov_params_val, data_object$randcov_Zs)

    # HW steps
    hwInv_val <- hwInv(smwInv_rand_val$SigInv, smwInv_rand_val$Sigldet, data_object$observed_index)

    SigInv <- hwInv_val$SigInv
    Sigldet <- hwInv_val$Sigldet
  } else {

    # making a covariance matrix
    cov_matrix_val_full <- cov_matrix(
      spcov_params_val, dist_matrix, randcov_params_val,
      data_object$randcov_Zs, data_object$partition_matrix, data_object$M
    )

    # subsetting for only observed
    cov_matrix_val <- cov_matrix_val_full[data_object$observed_index, data_object$observed_index, drop = FALSE]

    ## find cholesky
    chol_cov_matrix_val <- chol(forceSymmetric(cov_matrix_val))

    ## store value
    SigInv <- chol2inv(chol_cov_matrix_val)
    Sigldet <- 2 * sum(log(diag(chol_cov_matrix_val)))
  }

  # store value
  if (ldet) {
    return(list(SigInv = SigInv, Sigldet = Sigldet))
  } else {
    # need to clean up rest of code to not compute ldet
    return(list(SigInv = SigInv, Sigldet = NULL))
  }
}
