get_model_stats_spglm <- function(cov_est_object, data_object, estmethod) {



  # making a covariance matrix list
  cov_matrix_list <- get_cov_matrix_list(
    cov_est_object$spcov_params_val,
    cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val,
    data_object$randcov_list,
    data_object$partition_list
  )

  # find model components
  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)

  # cholesky products
  cholprods_list <- mapply(
    c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
    function(c, x, y) get_cholprods_glm(c, x, y),
    SIMPLIFY = FALSE
  )

  SigInv_list <- lapply(cholprods_list, function(x) x$SigInv)
  SigInv <- Matrix::bdiag(SigInv_list)
  SigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SigInv_X))

  # cov adjustment code
  # get cov beta hat (Xt Sig^-1 X)^-1 and beta hat (Xt Sig^-1 X)^-1 Xt Sig^-1 y
  invcov_betahat_list <- lapply(cholprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_X))

  invcov_betahat_sum <- Reduce("+", invcov_betahat_list)
  cov_betahat_noadjust <- chol2inv(chol(forceSymmetric(invcov_betahat_sum)))
  cov_betahat_noadjust_list <- rep(list(cov_betahat_noadjust), times = length(invcov_betahat_list))


  # find dispersion
  dispersion <- as.vector(cov_est_object$dispersion_params_val) # take class away

  # newton rhapson
  w_and_H <- get_w_and_H_spglm(data_object, dispersion,
                               SigInv_list, SigInv_X, cov_betahat_noadjust,
                               invcov_betahat_sum, estmethod,
                               ret_mHInv = TRUE)

  w <- w_and_H$w
  # H <- w_and_H$H

  # put w in cholprods
  w_list <- split(w, data_object$local_index)
  Xt_SigInv_w_list <- mapply(
    x = cholprods_list, w = w_list,
    function(x, w) crossprod(x$SqrtSigInv_X, forwardsolve(x$Sig_lowchol, w)),
    SIMPLIFY = FALSE
  )

  betahat_list <- mapply(
    l = cov_betahat_noadjust_list, r = Xt_SigInv_w_list,
    function(l, r) l %*% r,
    SIMPLIFY = FALSE
  )

  betahat <- as.numeric(cov_betahat_noadjust %*%
                          Reduce("+", Xt_SigInv_w_list))
  names(betahat) <- colnames(data_object$X_list[[1]])

  cov_betahat <- cov_betahat_adjust(
    invcov_betahat_list,
    betahat_list, betahat,
    cholprods_list, data_object,
    cov_est_object$spcov_params_val,
    cov_est_object$randcov_params_val,
    cov_betahat_noadjust, data_object$var_adjust
  )
  cov_betahat <- as.matrix(cov_betahat)
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  betawtsvarw <- wts_beta %*% w_and_H$mHInv %*% t(wts_beta)
  cov_betahat <- as.matrix(cov_betahat + betawtsvarw)

  rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat) <- colnames(data_object$X_list[[1]])

  # return coefficients
  coefficients <- get_coefficients_glm(betahat, cov_est_object$spcov_params_val,
                                       cov_est_object$dispersion_params_val, cov_est_object$randcov_params_val)

  # return fitted
  fitted <- get_fitted_glm(w, data_object)

}
