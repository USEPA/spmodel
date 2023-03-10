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
    eigenprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_eigenprods_glm_parallel)
    names(eigenprods_list) <- names(cov_matrix_list)
  } else {
    eigenprods_list <- mapply(
      c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, o = data_object$ones_list,
      function(c, x, y, o) get_eigenprods_glm(c, x, y, o),
      SIMPLIFY = FALSE
    )
  }

  SigInv_list <- lapply(eigenprods_list, function(x) x$SigInv)
  SigInv <- Matrix::bdiag(SigInv_list)
  SigInv_X <- do.call("rbind", lapply(eigenprods_list, function(x) x$SigInv_X))

  # cov adjustment code
  # get cov beta hat (Xt Sig^-1 X)^-1 and beta hat (Xt Sig^-1 X)^-1 Xt Sig^-1 y
  invcov_betahat_list <- lapply(eigenprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_X))

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

  # put w in eigenprods
  w_list <- split(w, sort(data_object$local_index))
  Xt_SigInv_w_list <- mapply(
    x = eigenprods_list, w = w_list,
    function(x, w) crossprod(x$SigInv_X, w),
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
    eigenprods_list, data_object,
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

  # return hat values
  hatvalues <- get_hatvalues_glm(w, X, data_object, dispersion)

  # return deviance i
  deviance_i <- get_deviance_glm(data_object$family, y, fitted$response, data_object$size, dispersion)
  deviance_i <- pmax(deviance_i, 0) # sometimes numerical instability can cause these to be slightly non-negative

  # storing relevant products
  SigInv_X_null <- do.call("rbind", lapply(eigenprods_list, function(x) x$SigInv_ones))
  ## lower chol %*% X
  SqrtSigInv_X_null <- do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_ones))
  # covariance of beta hat
  ## t(X) %*% sigma_inverse %*% X
  Xt_SigInv_X_null <- crossprod(SqrtSigInv_X_null, SqrtSigInv_X_null)
  ## t(X) %*% sigma_inverse %*% X)^(-1)
  Xt_SigInv_X_upchol_null <- chol(Xt_SigInv_X_null)
  cov_betahat_null <- chol2inv(Xt_SigInv_X_upchol_null)

  # newton rhapson
  w_and_H_null <- get_w_and_H_spglm(data_object, dispersion,
                                  SigInv_list, SigInv_X_null, cov_betahat_null, Xt_SigInv_X_null, estmethod)


  w_null <- w_and_H_null$w

  fitted_null <- get_fitted_glm(w_null, data_object)$response

  # return deviance i
  deviance_i_null <- get_deviance_glm(data_object$family, y, fitted_null, data_object$size, dispersion)
  deviance_i_null <- pmax(deviance_i_null, 0) # sometimes numerical instability can cause these to be slightly non-negative

  deviance <- sum(deviance_i)
  deviance_null <- sum(deviance_i_null)
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)

  # should always be non-negative
  pseudoR2 <- pmax(0, pseudoR2)
  # set null model R2 equal to zero (no covariates)
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }

  # return residuals
  residuals <- get_residuals_glm(w, y, data_object, deviance_i, hatvalues, dispersion)

  # return cooks distance
  cooks_distance <- get_cooks_distance_glm(residuals, hatvalues, data_object$p)

  # return variance covariance matrices
  vcov <- get_vcov_glm(cov_betahat) # note this is the adjusted one

  # reorder relevant quantities
  ## fitted
  fitted$response <- fitted$response[order(data_object$order)]
  names(fitted$response) <- data_object$observed_index
  fitted$link <- fitted$link[order(data_object$order)]
  names(fitted$link) <- data_object$observed_index
  hatvalues <- hatvalues[order(data_object$order)]
  names(hatvalues) <- data_object$observed_index
  residuals$response <- residuals$response[order(data_object$order)]
  names(residuals$response) <- data_object$observed_index
  residuals$deviance <- residuals$deviance[order(data_object$order)]
  names(residuals$deviance) <- data_object$observed_index
  residuals$pearson <- residuals$pearson[order(data_object$order)]
  names(residuals$pearson) <- data_object$observed_index
  residuals$standardized <- residuals$standardized[order(data_object$order)]
  names(residuals$standardized) <- data_object$observed_index
  cooks_distance <- cooks_distance[order(data_object$order)]
  names(cooks_distance) <- data_object$observed_index
  y <- y[order(data_object$order)]
  if (is.null(data_object$size)) {
    size <- NULL
  } else {
    size <- data_object$size[order(data_object$order)]
  }


  # return npar
  npar <- sum(unlist(lapply(cov_est_object$is_known, function(x) length(x) - sum(x)))) # could do sum(!x$is_known)

  # return list
  list(
    coefficients = coefficients,
    fitted = fitted,
    hatvalues = hatvalues,
    residuals = residuals,
    cooks_distance = cooks_distance,
    vcov = vcov,
    deviance = deviance,
    pseudoR2 = pseudoR2,
    npar = npar,
    w = w,
    y = y, # problems with model.response later
    size = size
  )
}

get_model_stats_spgautor <- function(cov_est_object, data_object, estmethod) {

  # cov_est_object$randcov_params_val is NULL if not added so won't affect downstream calculations
  # when random effects are not used

  cov_matrix_val <- cov_matrix(
    cov_est_object$spcov_params_val, cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val, data_object$randcov_Zs, data_object$partition_matrix, data_object$M
  )

  X <- data_object$X
  y <- data_object$y

  cov_matrix_obs_val <- cov_matrix_val[data_object$observed_index, data_object$observed_index, drop = FALSE]

  # getting cholesky products
  eigenprods <- get_eigenprods_glm(cov_matrix_obs_val, data_object$X, data_object$y, data_object$ones)

  SigInv <- eigenprods$SigInv
  # finding relevant quantities for likelihood
  SigInv_X <- SigInv %*% data_object$X
  Xt_SigInv_X <- crossprod(data_object$X, SigInv_X)
  Xt_SigInv_X_upchol <- chol(forceSymmetric(Xt_SigInv_X))
  cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

  # find dispersion
  dispersion <- as.vector(cov_est_object$dispersion_params_val) # take class away

  # newton rhapson
  w_and_H <- get_w_and_H_spgautor(data_object, dispersion,
                                  SigInv, SigInv_X, cov_betahat, Xt_SigInv_X, estmethod,
                                  ret_mHInv = TRUE)

  w <- w_and_H$w

  betahat <- as.numeric(tcrossprod(cov_betahat, SigInv_X) %*% w)
  names(betahat) <- colnames(data_object$X)

  cov_betahat <- as.matrix(cov_betahat)
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  betawtsvarw <- wts_beta %*% w_and_H$mHInv %*% t(wts_beta)
  cov_betahat <- as.matrix(cov_betahat + betawtsvarw)

  rownames(cov_betahat) <- colnames(data_object$X)
  colnames(cov_betahat) <- colnames(data_object$X)

  # return coefficients
  coefficients <- get_coefficients_glm(betahat, cov_est_object$spcov_params_val,
                                       cov_est_object$dispersion_params_val, cov_est_object$randcov_params_val)

  # return fitted
  fitted <- get_fitted_glm(w, data_object)

  # return hat values
  hatvalues <- get_hatvalues_glm(w, X, data_object, dispersion)

  # return deviance i
  deviance_i <- get_deviance_glm(data_object$family, y, fitted$response, data_object$size, dispersion)
  deviance_i <- pmax(deviance_i, 0) # sometimes numerical instability can cause these to be slightly non-negative

  # storing relevant products
  SigInv_X_null <- eigenprods_null$SigInv_X
  ## lower chol %*% X
  SqrtSigInv_X_null <- eigenprods_null$SqrtSigInv_X
  # covariance of beta hat
  ## t(X) %*% sigma_inverse %*% X
  Xt_SigInv_X_null <- crossprod(SqrtSigInv_X_null, SqrtSigInv_X_null)
  ## t(X) %*% sigma_inverse %*% X)^(-1)
  Xt_SigInv_X_upchol_null <- chol(Xt_SigInv_X_null)
  cov_betahat_null <- chol2inv(Xt_SigInv_X_upchol_null)

  # newton rhapson
  w_and_H_null <- get_w_and_H_spgautor(data_object, dispersion,
                                    SigInv_null, SigInv_X_null, cov_betahat_null, Xt_SigInv_X_null, estmethod)


  w_null <- w_and_H_null$w

  fitted_null <- get_fitted_glm(w_null, data_object)$response

  # return deviance i
  deviance_i_null <- get_deviance_glm(data_object$family, y, fitted_null, data_object$size, dispersion)
  deviance_i_null <- pmax(deviance_i_null, 0) # sometimes numerical instability can cause these to be slightly non-negative

  deviance <- sum(deviance_i)
  deviance_null <- sum(deviance_i_null)
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)

  # should always be non-negative
  pseudoR2 <- pmax(0, pseudoR2)
  # set null model R2 equal to zero (no covariates)
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }

  # return residuals
  residuals <- get_residuals_glm(w, y, data_object, deviance_i, hatvalues, dispersion)

  # return cooks distance
  cooks_distance <- get_cooks_distance_glm(residuals, hatvalues, data_object$p)

  # give names
  names(fitted$response) <- data_object$observed_index
  names(fitted$link) <- data_object$observed_index
  names(hatvalues) <- data_object$observed_index
  names(residuals$response) <- data_object$observed_index
  names(residuals$deviance) <- data_object$observed_index
  names(residuals$pearson) <- data_object$observed_index
  names(residuals$standardized) <- data_object$observed_index
  names(cooks_distance) <- data_object$observed_index

  # return variance covariance matrices
  vcov <- get_vcov_glm(cov_betahat) # note this is the adjusted one

  # return npar
  npar <- sum(unlist(lapply(cov_est_object$is_known, function(x) length(x) - sum(x)))) # could do sum(!x$is_known)

  # return list
  list(
    coefficients = coefficients,
    fitted = fitted,
    hatvalues = hatvalues,
    residuals = residuals,
    cooks_distance = cooks_distance,
    vcov = vcov,
    deviance = deviance,
    pseudoR2 = pseudoR2,
    npar = npar,
    w = w,
    y = y,
    size = data_object$size
  )

}
