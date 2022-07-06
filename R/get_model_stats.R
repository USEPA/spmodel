get_model_stats_splm <- function(cov_est_object, data_object, estmethod) {


  # making a covariance matrix list
  cov_matrix_list <- get_cov_matrix_list(
    cov_est_object$spcov_params_val,
    cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val,
    data_object$randcov_list,
    data_object$partition_list
  )



  # cholesky products
  cholprods_list <- mapply(
    c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
    function(c, x, y) get_cholprods(c, x, y),
    SIMPLIFY = FALSE
  )


  # get cov beta hat (Xt Sig^-1 X)^-1 and beta hat (Xt Sig^-1 X)^-1 Xt Sig^-1 y
  invcov_betahat_list <- lapply(cholprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_X))

  invcov_betahat_sum <- Reduce("+", invcov_betahat_list)
  cov_betahat_noadjust <- chol2inv(chol(forceSymmetric(invcov_betahat_sum)))
  cov_betahat_noadjust_list <- rep(list(cov_betahat_noadjust), times = length(invcov_betahat_list))
  # cov_betahat_noadjust_list <- lapply(invcov_betahat_list, function(x) chol2inv(chol(forceSymmetric(x))))
  # not PD / wrong answer when partition factor in fixed effects (no variation in fixed effects)
  # need to figure out what to do with pooled cov adjust

  Xt_SigInv_y_list <- lapply(cholprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_y))

  betahat_list <- mapply(
    l = cov_betahat_noadjust_list, r = Xt_SigInv_y_list,
    function(l, r) l %*% r,
    SIMPLIFY = FALSE
  )

  betahat <- as.numeric(cov_betahat_noadjust %*%
    Reduce("+", Xt_SigInv_y_list))
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
  rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat) <- colnames(data_object$X_list[[1]])


  # return coefficients
  coefficients <- get_coefficients(betahat, cov_est_object$spcov_params_val, cov_est_object$randcov_params_val)

  # return fitted
  fitted <- get_fitted_splm(
    betahat, cov_est_object$spcov_params_val, data_object,
    cholprods_list, cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val
  )

  # return hat values
  hatvalues <- as.numeric(unlist(lapply(cholprods_list, function(x) get_hatvalues(cov_betahat, x$SqrtSigInv_X))))

  # return residuals
  residuals <- get_residuals_splm(betahat, data_object, cholprods_list, hatvalues)

  # return cooks distance
  cooks_distance <- get_cooks_distance(residuals, hatvalues, data_object$p)

  # reorder relevant quantities
  ## fitted
  fitted$response <- fitted$response[order(data_object$order)]
  fitted$spcov$de <- fitted$spcov$de[order(data_object$order)]
  fitted$spcov$ie <- fitted$spcov$ie[order(data_object$order)]
  hatvalues <- hatvalues[order(data_object$order)]
  residuals$raw <- residuals$raw[order(data_object$order)]
  residuals$pearson <- residuals$pearson[order(data_object$order)]
  residuals$standardized <- residuals$standardized[order(data_object$order)]
  cooks_distance <- cooks_distance[order(data_object$order)]

  # return variance covariance matrices
  vcov <- get_vcov(cov_betahat)

  # return deviance
  deviance <- as.numeric(crossprod(residuals$pearson, residuals$pearson))

  # generalized r squared
  ## generalized r squared (1 - deviance full / deviance reduced (mean only))
  ## for normal data is 1 - (y - x beta)t Sigma inv (y - x beta) / (y - muhat)t Sigma inv (y - muhat)
  ## where muhat = (1t Sigma inv 1) inv 1t Sigma inv y
  ## create ones vector
  ## muhat
  SqrtSigInv_ones <- as.numeric(do.call(
    "rbind",
    mapply(
      c = cholprods_list, o = data_object$ones_list,
      function(c, o) forwardsolve(c$Sig_lowchol, o), SIMPLIFY = FALSE
    )
  ))
  cov_muhat <- 1 / crossprod(SqrtSigInv_ones, SqrtSigInv_ones)
  SqrtSigInv_y <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_y))
  muhat <- cov_muhat * crossprod(SqrtSigInv_ones, SqrtSigInv_y)
  SqrtSigInv_rmuhat <- as.numeric(do.call(
    "rbind",
    mapply(
      c = cholprods_list, y = data_object$y_list,
      function(c, y) forwardsolve(c$Sig_lowchol, y - rep(muhat, length(y))), SIMPLIFY = FALSE
    )
  ))
  ### reduced model
  deviance_null <- as.numeric(crossprod(SqrtSigInv_rmuhat, SqrtSigInv_rmuhat))
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)

  # set null model R2 equal to zero (no covariates)
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }


  ## empirical semivariogram
  if (estmethod == "sv-wls") {
    esv_val <- cov_est_object$esv
  } else {
    esv_val <- NULL
  }


  # npar
  if (estmethod == "ml") {
    p_theta_spcov <- length(cov_est_object$is_known$spcov) - sum(cov_est_object$is_known$spcov)
    p_theta_randcov <- length(cov_est_object$is_known$randcov) - sum(cov_est_object$is_known$randcov)
    p_theta <- p_theta_spcov + p_theta_randcov
    p_beta <- length(coefficients$fixed)
    npar <- p_theta + p_beta
  } else {
    p_theta_spcov <- length(cov_est_object$is_known$spcov) - sum(cov_est_object$is_known$spcov)
    p_theta_randcov <- length(cov_est_object$is_known$randcov) - sum(cov_est_object$is_known$randcov)
    p_theta <- p_theta_spcov + p_theta_randcov
    p_beta <- 0
    npar <- p_theta + p_beta
  }

  list(
    coefficients = coefficients,
    fitted = fitted,
    hatvalues = hatvalues,
    residuals = residuals,
    cooks_distance = cooks_distance,
    vcov = vcov,
    deviance = deviance,
    pseudoR2 = pseudoR2,
    npar = npar
  )
}


get_model_stats_spautor <- function(cov_est_object, data_object, estmethod) {
  # cov_est_object$randcov_params_val is NULL if not added so won't affect downstream calculations
  # when random effects are not used

  # USE GLOGLIK_PRODUCT CODE TO AVOID CHOLESKY IF ABLE (IE = 0, NO RANDOM OR PARTITIONING)
  # REQUIRES REWRITING DIAGNOSTICS
  cov_matrix_val <- cov_matrix(
    cov_est_object$spcov_params_val, cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val, data_object$randcov_Zs, data_object$partition_matrix, data_object$M
  )

  cov_matrix_obs_val <- cov_matrix_val[data_object$observed_index, data_object$observed_index, drop = FALSE]

  # getting cholesky products
  cholprods <- get_cholprods(cov_matrix_obs_val, data_object$X, data_object$y)

  # get cov beta hat (Xt Sig^-1 X)^-1
  cov_betahat <- as.matrix(chol2inv(chol(forceSymmetric(crossprod(cholprods$SqrtSigInv_X, cholprods$SqrtSigInv_X)))))
  rownames(cov_betahat) <- colnames(data_object$X)
  colnames(cov_betahat) <- colnames(data_object$X)

  # get betahat (Xt Sig^-1 X)^-1 Xt Sig^-1 y
  betahat <- cov_betahat %*% crossprod(cholprods$SqrtSigInv_X, cholprods$SqrtSigInv_y)
  betahat <- as.numeric(betahat)
  names(betahat) <- colnames(data_object$X)

  # return coefficients
  coefficients <- get_coefficients(betahat, cov_est_object$spcov_params_val, cov_est_object$randcov_params_val)

  # return fitted
  fitted <- get_fitted_spautor(
    betahat, cov_est_object$spcov_params_val, data_object, cholprods,
    cov_est_object$randcov_params_val
  )

  # return hat values
  hatvalues <- get_hatvalues(cov_betahat, cholprods$SqrtSigInv_X)

  # return residuals
  residuals <- get_residuals_spautor(betahat, data_object$X, data_object$y, cholprods, hatvalues)

  # return cooks distance
  cooks_distance <- get_cooks_distance(residuals, hatvalues, data_object$p)

  # return variance covariance matrices
  vcov <- get_vcov(cov_betahat)

  # return deviance
  deviance <- as.numeric(crossprod(residuals$pearson, residuals$pearson))

  # generalized r squared
  ## generalized r squared (1 - deviance full / deviance reduced (mean only))
  ## for normal data is 1 - (y - x beta)t Sigma inv (y - x beta) / (y - muhat)t Sigma inv (y - muhat)
  ## where muhat = (1t Sigma inv 1) inv 1t Sigma inv y
  ## create ones vector
  ## muhat
  SqrtSigInv_ones <- forwardsolve(cholprods$Sig_lowchol, data_object$ones)
  cov_muhat <- 1 / crossprod(SqrtSigInv_ones, SqrtSigInv_ones)
  muhat <- cov_muhat * crossprod(SqrtSigInv_ones, cholprods$SqrtSigInv_y)
  SqrtSigInv_rmuhat <- forwardsolve(cholprods$Sig_lowchol, data_object$y - rep(muhat, data_object$n))
  ### reduced model
  deviance_null <- as.numeric(crossprod(SqrtSigInv_rmuhat, SqrtSigInv_rmuhat))
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)

  # set null model R2 equal to zero (no covariates)
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }

  # npar
  if (estmethod == "ml") {
    p_theta_spcov <- length(cov_est_object$is_known$spcov) - sum(cov_est_object$is_known$spcov)
    p_theta_randcov <- length(cov_est_object$is_known$randcov) - sum(cov_est_object$is_known$randcov)
    p_theta <- p_theta_spcov + p_theta_randcov
    p_beta <- length(coefficients$fixed)
    npar <- p_theta + p_beta
  } else {
    p_theta_spcov <- length(cov_est_object$is_known$spcov) - sum(cov_est_object$is_known$spcov)
    p_theta_randcov <- length(cov_est_object$is_known$randcov) - sum(cov_est_object$is_known$randcov)
    p_theta <- p_theta_spcov + p_theta_randcov
    p_beta <- 0
    npar <- p_theta + p_beta
  }

  list(
    coefficients = coefficients,
    fitted = fitted,
    hatvalues = hatvalues,
    residuals = residuals,
    cooks_distance = cooks_distance,
    vcov = vcov,
    deviance = deviance,
    pseudoR2 = pseudoR2,
    npar = npar
  )
}
