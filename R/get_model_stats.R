get_model_stats_splm <- function(cov_est_object, data_object, estmethod) {


  # making a covariance matrix list
  cov_matrix_list <- get_cov_matrix_list(
    cov_est_object$spcov_params_val,
    cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val,
    data_object$randcov_list,
    data_object$partition_list,
    diagtol = data_object$diagtol
  )



  # eigen products
  if (data_object$parallel) {
    cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
      cluster_list_element <- list(
        c = cov_matrix_list[[l]],
        x = data_object$X_list[[l]],
        y = data_object$y_list[[l]],
        o = data_object$ones_list[[l]]
      )
    })
    eigenprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_eigenprods_parallel)
    names(eigenprods_list) <- names(cov_matrix_list)
  } else {
    eigenprods_list <- mapply(
      c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, o = data_object$ones_list,
      function(c, x, y, o) get_eigenprods(c, x, y, o),
      SIMPLIFY = FALSE
    )
  }

  # get cov beta hat (Xt Sig^-1 X)^-1 and beta hat (Xt Sig^-1 X)^-1 Xt Sig^-1 y
  invcov_betahat_list <- lapply(eigenprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_X))

  invcov_betahat_sum <- Reduce("+", invcov_betahat_list)
  cov_betahat_noadjust <- chol2inv(chol(forceSymmetric(invcov_betahat_sum)))
  cov_betahat_noadjust_list <- rep(list(cov_betahat_noadjust), times = length(invcov_betahat_list))

  Xt_SigInv_y_list <- lapply(eigenprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_y))

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
    eigenprods_list, data_object,
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
    eigenprods_list, cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val
  )

  # return hat values
  hatvalues <- as.numeric(unlist(lapply(eigenprods_list, function(x) get_hatvalues(cov_betahat_noadjust, x$SqrtSigInv_X))))

  # return residuals
  residuals <- get_residuals_splm(betahat, data_object, eigenprods_list, hatvalues)

  # return cooks distance
  cooks_distance <- get_cooks_distance(residuals, hatvalues, data_object$p)

  # reorder relevant quantities
  ## fitted
  fitted$response <- fitted$response[order(data_object$order)]
  names(fitted$response) <- data_object$observed_index
  fitted$spcov$de <- fitted$spcov$de[order(data_object$order)]
  names(fitted$spcov$de) <- data_object$observed_index
  fitted$spcov$ie <- fitted$spcov$ie[order(data_object$order)]
  names(fitted$spcov$ie) <- data_object$observed_index
  hatvalues <- hatvalues[order(data_object$order)]
  names(hatvalues) <- data_object$observed_index
  residuals$response <- residuals$response[order(data_object$order)]
  names(residuals$response) <- data_object$observed_index
  residuals$pearson <- residuals$pearson[order(data_object$order)]
  names(residuals$pearson) <- data_object$observed_index
  residuals$standardized <- residuals$standardized[order(data_object$order)]
  names(residuals$standardized) <- data_object$observed_index
  cooks_distance <- cooks_distance[order(data_object$order)]
  names(cooks_distance) <- data_object$observed_index

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
  SqrtSigInv_ones <- as.numeric(do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_ones)))
  cov_muhat <- 1 / crossprod(SqrtSigInv_ones, SqrtSigInv_ones)
  SqrtSigInv_y <- do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_y))
  muhat <- as.vector(cov_muhat * crossprod(SqrtSigInv_ones, SqrtSigInv_y))
  SqrtSigInv_rmuhat <- as.numeric(do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_y - x$SqrtSigInv_ones * muhat)))
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
  p_theta_spcov <- length(cov_est_object$is_known$spcov) - sum(cov_est_object$is_known$spcov)
  p_theta_randcov <- length(cov_est_object$is_known$randcov) - sum(cov_est_object$is_known$randcov)
  npar <- p_theta_spcov + p_theta_randcov

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

get_model_stats_splm_iid <- function(cov_est_object, data_object, estmethod) {


  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)
  qr_val <- qr(X)
  R_val <- qr.R(qr_val)
  s2 <- cov_est_object$spcov_params_val[["ie"]]
  cor_betahat <- chol2inv(chol(crossprod(R_val, R_val)))
  cov_betahat <- s2 * cor_betahat
  betahat <- as.numeric(backsolve(R_val, qr.qty(qr_val, y)))
  names(betahat) <- colnames(data_object$X_list[[1]])
  cov_betahat <- as.matrix(cov_betahat)
  rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat) <- colnames(data_object$X_list[[1]])
  fitted <- X %*% betahat
  resids <- y - fitted


  # return coefficients
  coefficients <- get_coefficients(betahat, cov_est_object$spcov_params_val, cov_est_object$randcov_params_val)

  # return fitted
  fitted <- list(
    response = as.numeric(fitted),
    spcov = list(de = as.numeric(rep(0, length(y))), ie = as.numeric(resids)),
    randcov = NULL
  )


  # return hat values
  hatvalues <- diag(X %*% tcrossprod(cor_betahat, X))
  # return residuals
  residuals <- list(
    response = as.numeric(resids),
    pearson = 1 / sqrt(s2) * resids
  )
  residuals$standardized <- residuals$pearson / sqrt(1 - hatvalues)

  # return cooks distance
  cooks_distance <- get_cooks_distance(residuals, hatvalues, data_object$p)

  # reorder relevant quantities
  ## fitted
  fitted$response <- fitted$response[order(data_object$order)]
  names(fitted$response) <- data_object$observed_index
  fitted$spcov$de <- fitted$spcov$de[order(data_object$order)]
  names(fitted$spcov$de) <- data_object$observed_index
  fitted$spcov$ie <- fitted$spcov$ie[order(data_object$order)]
  names(fitted$spcov$ie) <- data_object$observed_index
  hatvalues <- hatvalues[order(data_object$order)]
  names(hatvalues) <- data_object$observed_index
  residuals$response <- residuals$response[order(data_object$order)]
  names(residuals$response) <- data_object$observed_index
  residuals$pearson <- residuals$pearson[order(data_object$order)]
  names(residuals$pearson) <- data_object$observed_index
  residuals$standardized <- residuals$standardized[order(data_object$order)]
  names(residuals$standardized) <- data_object$observed_index
  cooks_distance <- cooks_distance[order(data_object$order)]
  names(cooks_distance) <- data_object$observed_index

  # return variance covariance matrices
  vcov <- get_vcov(cov_betahat)

  # return deviance
  deviance <- as.numeric(crossprod(residuals$pearson, residuals$pearson))
  muhat <- mean(y)
  pearson_null <- 1 / sqrt(s2) * (y - muhat)
  deviance_null <- as.numeric(crossprod(pearson_null, pearson_null))
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
  p_theta_spcov <- length(cov_est_object$is_known$spcov) - sum(cov_est_object$is_known$spcov)
  p_theta_randcov <- length(cov_est_object$is_known$randcov) - sum(cov_est_object$is_known$randcov)
  npar <- p_theta_spcov + p_theta_randcov

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

  cov_matrix_val <- cov_matrix(
    cov_est_object$spcov_params_val, cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val, data_object$randcov_Zs, data_object$partition_matrix, data_object$M
  )

  cov_matrix_obs_val <- cov_matrix_val[data_object$observed_index, data_object$observed_index, drop = FALSE]

  # getting cholesky products
  eigenprods <- get_eigenprods(cov_matrix_obs_val, data_object$X, data_object$y, data_object$ones)

  # get cov beta hat (Xt Sig^-1 X)^-1
  cov_betahat <- as.matrix(chol2inv(chol(forceSymmetric(crossprod(eigenprods$SqrtSigInv_X, eigenprods$SqrtSigInv_X)))))
  rownames(cov_betahat) <- colnames(data_object$X)
  colnames(cov_betahat) <- colnames(data_object$X)

  # get betahat (Xt Sig^-1 X)^-1 Xt Sig^-1 y
  betahat <- cov_betahat %*% crossprod(eigenprods$SqrtSigInv_X, eigenprods$SqrtSigInv_y)
  betahat <- as.numeric(betahat)
  names(betahat) <- colnames(data_object$X)

  # return coefficients
  coefficients <- get_coefficients(betahat, cov_est_object$spcov_params_val, cov_est_object$randcov_params_val)

  # return fitted
  fitted <- get_fitted_spautor(
    betahat, cov_est_object$spcov_params_val, data_object, eigenprods,
    cov_est_object$randcov_params_val
  )

  # return hat values
  hatvalues <- get_hatvalues(cov_betahat, eigenprods$SqrtSigInv_X)

  # return residuals
  residuals <- get_residuals_spautor(betahat, data_object$X, data_object$y, eigenprods, hatvalues)

  # return cooks distance
  cooks_distance <- get_cooks_distance(residuals, hatvalues, data_object$p)

  # give names
  names(fitted$response) <- data_object$observed_index
  names(fitted$spcov$de) <- data_object$observed_index
  names(fitted$spcov$ie) <- data_object$observed_index
  names(hatvalues) <- data_object$observed_index
  names(residuals$response) <- data_object$observed_index
  names(residuals$pearson) <- data_object$observed_index
  names(residuals$standardized) <- data_object$observed_index
  names(cooks_distance) <- data_object$observed_index

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
  SqrtSigInv_ones <- eigenprods$SqrtSigInv_ones
  cov_muhat <- 1 / crossprod(SqrtSigInv_ones, SqrtSigInv_ones)
  muhat <- as.vector(cov_muhat * crossprod(SqrtSigInv_ones, eigenprods$SqrtSigInv_y))
  SqrtSigInv_rmuhat <- eigenprods$SqrtSigInv_y - eigenprods$SqrtSigInv_ones * muhat
  ### reduced model
  deviance_null <- as.numeric(crossprod(SqrtSigInv_rmuhat, SqrtSigInv_rmuhat))
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)

  # set null model R2 equal to zero (no covariates)
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }

  # npar
  p_theta_spcov <- length(cov_est_object$is_known$spcov) - sum(cov_est_object$is_known$spcov)
  p_theta_randcov <- length(cov_est_object$is_known$randcov) - sum(cov_est_object$is_known$randcov)
  npar <- p_theta_spcov + p_theta_randcov

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
