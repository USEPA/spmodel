#' Assemble fitted-model statistics for an \code{spglm()} model
#'
#' @param cov_est_object A fitted covariance estimation object
#' @param data_object The data object
#' @param estmethod The estimation method
#'
#' @return A list with the fixed effect coefficients, fitted values, hat
#'   values, deviance residuals, Cook's distance, covariance matrix of the
#'   fixed effects, deviance, pseudo R-squared, number of estimated covariance
#'   parameters, latent predictor vector \code{w}, response vector, and
#'   binomial size vector
#'
#' @noRd
get_model_stats_spglm <- function(cov_est_object, data_object, estmethod) {
  # GLM-type models don't have a closed-form likelihood the way splm() does,
  # so this function alternates: (1) treat the current betahat as fixed and
  # find the latent link-scale predictor w that maximizes the Laplace
  # approximation via Newton-Raphson (get_w_and_H_spglm()), then (2) treat w
  # as fixed and re-estimate betahat by (weighted) least squares -- repeated
  # once here for the full model and once more below for the intercept-only
  # ("null") model needed for deviance/pseudo R^2
  # making a covariance matrix list
  cov_matrix_list <- get_cov_matrix_list(
    cov_est_object$spcov_params_val,
    cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val,
    data_object$randcov_list,
    data_object$partition_list,
    diagtol = data_object$diagtol
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
  # find w (the latent link-scale predictor) that maximizes the Laplace
  # approximation to the likelihood, along with mHInv, the inverse negative
  # Hessian at that maximum -- mHInv quantifies the extra uncertainty from
  # estimating w rather than observing it directly, used below in betawtsvarw
  w_and_H <- get_w_and_H_spglm(data_object, dispersion,
    SigInv_list, SigInv_X, cov_betahat_noadjust,
    invcov_betahat_sum, estmethod,
    ret_mHInv = TRUE
  )

  w <- as.vector(w_and_H$w)
  # H <- w_and_H$H

  # put w in eigenprods
  # local_index tracks which partition each observation belongs to; sorting
  # by it and splitting recovers the same per-partition grouping used to
  # build eigenprods_list, so w can be combined partition-by-partition below
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
  # additional variance in betahat induced by w itself being estimated
  # (rather than known), propagated through via the delta method using mHInv
  betawtsvarw <- wts_beta %*% w_and_H$mHInv %*% t(wts_beta)

  cov_betahat_uncorrected <- cov_betahat # save uncorrected cov beta hat
  rownames(cov_betahat_uncorrected) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat_uncorrected) <- colnames(data_object$X_list[[1]])

  cov_betahat <- as.matrix(cov_betahat + betawtsvarw)
  rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat) <- colnames(data_object$X_list[[1]])

  # return coefficients
  coefficients <- get_coefficients_glm(
    betahat, cov_est_object$spcov_params_val,
    cov_est_object$dispersion_params_val, cov_est_object$randcov_params_val
  )

  # return fitted
  fitted <- get_fitted_spglm(
    w_list, betahat, cov_est_object$spcov_params_val, data_object,
    eigenprods_list, cov_est_object$dist_matrix_list,
    cov_est_object$randcov_params_val
  )

  # return hat values
  hatvalues <- get_hatvalues_glm(w, X, data_object, dispersion)

  # return deviance i
  deviance_i <- get_deviance_glm(data_object$family, y, fitted$response, data_object$size, dispersion)
  deviance_i <- pmax(deviance_i, 0) # sometimes numerical instability can cause these to be slightly negative

  # below refits an intercept-only ("null") model the same way as above (its
  # own Newton-Raphson for w, its own betahat) so deviance_i_null can be
  # compared against deviance_i to form pseudo R^2 further down
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
  w_and_H_null <- get_w_and_H_spglm(
    data_object, dispersion,
    SigInv_list, SigInv_X_null, cov_betahat_null, Xt_SigInv_X_null, estmethod
  )


  w_null <- as.vector(w_and_H_null$w)

  fitted_null <- get_fitted_null(w_null, data_object)

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
  vcov <- get_vcov_glm(cov_betahat, cov_betahat_uncorrected) # note first argument is the adjusted one

  # local estimation partitions the data, so these vectors are currently
  # ordered by partition rather than original row order; data_object$order
  # maps partition order back to the original data order
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
  # counts only covariance/dispersion parameters actually estimated, i.e.
  # excludes any the user fixed at a known value (is_known == TRUE)
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

#' Assemble fitted-model statistics for an \code{spgautor()} model
#'
#' @param cov_est_object A fitted covariance estimation object
#' @param data_object The data object
#' @param estmethod The estimation method
#'
#' @return The same statistics as \code{get_model_stats_spglm()}, computed from the
#'   full neighborhood covariance matrix (\code{cov_est_object$randcov_params_val}
#'   is \code{NULL} when random effects are not used, so it does not affect
#'   downstream calculations in that case)
#'
#' @noRd
get_model_stats_spgautor <- function(cov_est_object, data_object, estmethod) {
  # spgautor() models a fixed neighborhood structure (no big-data local
  # partitioning), so this mirrors get_model_stats_spglm() but works with a
  # single covariance matrix/eigenprods object rather than per-partition lists
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
  # find w (the latent link-scale predictor) maximizing the Laplace
  # approximation, and mHInv (inverse negative Hessian at that maximum),
  # which quantifies the extra uncertainty from estimating w -- used below
  # in betawtsvarw
  w_and_H <- get_w_and_H_spgautor(data_object, dispersion,
    SigInv, SigInv_X, cov_betahat, Xt_SigInv_X, estmethod,
    ret_mHInv = TRUE
  )

  w <- as.vector(w_and_H$w) # remember this is w - offset so must add offset in relevant places later to match glm

  betahat <- as.numeric(tcrossprod(cov_betahat, SigInv_X) %*% w)
  names(betahat) <- colnames(data_object$X)

  cov_betahat <- as.matrix(cov_betahat)
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  # additional variance in betahat induced by w being estimated rather than
  # known, propagated via the delta method using mHInv
  betawtsvarw <- wts_beta %*% w_and_H$mHInv %*% t(wts_beta)

  cov_betahat_uncorrected <- cov_betahat # save uncorrected cov beta hat
  rownames(cov_betahat_uncorrected) <- colnames(data_object$X)
  colnames(cov_betahat_uncorrected) <- colnames(data_object$X)

  cov_betahat <- as.matrix(cov_betahat + betawtsvarw)
  rownames(cov_betahat) <- colnames(data_object$X)
  colnames(cov_betahat) <- colnames(data_object$X)

  # return coefficients
  coefficients <- get_coefficients_glm(
    betahat, cov_est_object$spcov_params_val,
    cov_est_object$dispersion_params_val, cov_est_object$randcov_params_val
  )

  # return fitted
  fitted <- get_fitted_spgautor(
    w, betahat, cov_est_object$spcov_params_val, data_object, eigenprods,
    randcov_params = cov_est_object$randcov_params_val
  )

  # return hat values
  hatvalues <- get_hatvalues_glm(w, X, data_object, dispersion)

  # return deviance i
  deviance_i <- get_deviance_glm(data_object$family, y, fitted$response, data_object$size, dispersion)
  deviance_i <- pmax(deviance_i, 0) # sometimes numerical instability can cause these to be slightly non-negative

  # below refits an intercept-only ("null") model the same way as above so
  # deviance_i_null can be compared against deviance_i to form pseudo R^2
  # storing relevant products
  SigInv_X_null <- eigenprods$SigInv_ones
  ## lower chol %*% X
  SqrtSigInv_X_null <- eigenprods$SqrtSigInv_ones
  # covariance of beta hat
  ## t(X) %*% sigma_inverse %*% X
  Xt_SigInv_X_null <- crossprod(SqrtSigInv_X_null, SqrtSigInv_X_null)
  ## t(X) %*% sigma_inverse %*% X)^(-1)
  Xt_SigInv_X_upchol_null <- chol(Xt_SigInv_X_null)
  cov_betahat_null <- chol2inv(Xt_SigInv_X_upchol_null)

  # newton rhapson
  w_and_H_null <- get_w_and_H_spgautor(
    data_object, dispersion,
    SigInv, SigInv_X_null, cov_betahat_null, Xt_SigInv_X_null, estmethod
  )


  w_null <- as.vector(w_and_H_null$w)

  fitted_null <- get_fitted_null(w_null, data_object)

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
  vcov <- get_vcov_glm(cov_betahat, cov_betahat_uncorrected) # note first argument is the adjusted one

  # return npar
  # counts only covariance/dispersion parameters actually estimated, i.e.
  # excludes any the user fixed at a known value (is_known == TRUE)
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
