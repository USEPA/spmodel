#' @param type The scale (\code{response} or \code{link}) of predictions obtained
#'   when \code{cv_predict = TRUE} and using \code{spglm()} or \code{spgautor} objects.
#' @param delta A logical indicating whether to return delta method standard errors
#' on the response scale when \code{se.fit = TRUE} and \code{type = "response"}. The default is \code{FALSE}.
#' @rdname loocv
#' @method loocv spglm
#' @order 4
#' @export
loocv.spglm <- function(object, cv_predict = FALSE, type = c("link", "response"), se.fit = FALSE, delta = FALSE, local, ...) {
  # match type argument so the two display
  type <- match.arg(type)

  if (missing(local)) local <- NULL

  # local prediction list

  # local stuff
  if (is.null(local)) {
    if (object$n > 5000) {
      local <- TRUE
      message("Because the sample size exceeds 5000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun loocv() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }
  local_list <- get_local_list_prediction(local)

  y <- object$y

  if (local_list$method == "all") {
    # exact LOO: build the covariance matrix (on the latent link scale) and
    # its Cholesky-based inverse once, then reuse for every held-out row via
    # get_loocv_glm()'s closed-form update instead of refitting n times
    cov_matrix_val <- covmatrix(object)
    X <- model.matrix(object)
    cholprods <- get_cholprods_glm(cov_matrix_val, X, y)
    # actually need inverse because of HW blocking
    SigInv <- chol2inv(cholprods$Sig_lowchol)
    SigInv_X <- backsolve(t(cholprods$Sig_lowchol), cholprods$SqrtSigInv_X)

    # find products
    Xt_SigInv_X <- crossprod(X, SigInv_X)
    Xt_SigInv_X_upchol <- base::chol(Xt_SigInv_X) # or Matrix::forceSymmetric()
    cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

    # glm stuff
    dispersion <- as.vector(coef(object, type = "dispersion")) # take class away
    w_linpred <- fitted(object, type = "link") # offset included
    size <- object$size

    # The kriging update below covaries via Sigma, which describes the
    # offset-free latent process, so w is stripped of the offset here and each
    # held-out row's own offset is added back to its prediction afterwards.
    # get_d()/get_D() differentiate the data model and stay on w_linpred;
    # see w_offset_free()
    model_offset <- model.offset(model.frame(object))
    w <- w_offset_free(w_linpred, model_offset)

    # some products
    SigInv_w <- SigInv %*% w
    wX <- cbind(w, X)
    SigInv_wX <- cbind(SigInv_w, SigInv_X)

    # find H stuff: D and Ptheta are the two pieces of the (negative) Hessian
    # of the penalized quasi-likelihood used to fit the latent link-scale
    # values w, so H below approximates the curvature needed to propagate
    # estimation uncertainty in w into the leave-one-out updates
    wts_beta <- tcrossprod(cov_betahat, SigInv_X)
    Ptheta <- SigInv - SigInv_X %*% wts_beta
    d <- get_d(object$family, w_linpred, y, size, dispersion)
    # and then the gradient vector
    # g <-  d - Ptheta %*% w
    # Next, compute H
    D <- get_D(object$family, w_linpred, y, size, dispersion)
    H <- D - Ptheta
    mHinv <- solve(-H) # chol2inv(chol(Matrix::forceSymmetric(-H))) # solve(-H)

    # parallel stuff
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv_glm,
        Sig = cov_matrix_val,
        SigInv = SigInv, Xmat = X, w = w, wX = wX,
        SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
      )
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), get_loocv_glm,
        Sig = cov_matrix_val,
        SigInv = SigInv, Xmat = X, w = as.matrix(w, ncol = 1), wX = wX,
        SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
      )
    }
    cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    # get_loocv_glm() predicted the offset-free latent process, so give each
    # held-out row its own offset back (the standard errors are unaffected,
    # since the offset is a known constant shift)
    if (!is.null(model_offset)) {
      cv_predict_val <- cv_predict_val + as.vector(model_offset)
    }
    if (se.fit) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  } else {
    # precompute everything that is invariant across held-out rows (see the
    # analogous comment in loocv.splm()) so the per-row worker
    # (loocv_local_glm()) can slice arrays directly instead of calling
    # predict(), which previously re-derived every one of these from
    # object$formula/object$obdata on every single row
    spcov_params_val <- coef(object, type = "spcov")
    dispersion_params_val <- as.vector(coef(object, type = "dispersion"))
    randcov_params_val <- coef(object, type = "randcov")

    xcoord <- object$xcoord
    ycoord <- object$ycoord
    obdata_pred <- object$obdata
    if (object$dim_coords == 1) {
      obdata_pred[[ycoord]] <- 0
    }
    if (object$anisotropy) {
      aniscoords <- transform_anis(obdata_pred, xcoord, ycoord,
        rotate = spcov_params_val[["rotate"]], scale = spcov_params_val[["scale"]]
      )
      obdata_pred[[xcoord]] <- aniscoords$xcoord_val
      obdata_pred[[ycoord]] <- aniscoords$ycoord_val
    }

    Xmat_full <- model.matrix(object)
    w_full <- fitted(object, type = "link")
    model_offset_full <- model.offset(model.frame(object))
    size_full <- object$size
    betahat <- coefficients(object)
    cov_betahat <- vcov(object, var_correct = FALSE)

    extra_randcov_list <- get_extra_randcov_list(object, object$obdata, newdata = object$obdata)
    extra_partition_list <- get_extra_partition_list(object, object$obdata, newdata = object$obdata)

    loocv_context <- list(
      obdata_pred = obdata_pred, Xmat_full = Xmat_full, y_full = y, w_full = w_full,
      model_offset_full = model_offset_full, size_full = size_full,
      xcoord = xcoord, ycoord = ycoord, spcov_params_val = spcov_params_val, random = object$random,
      randcov_params_val = randcov_params_val, partition_factor = object$partition_factor,
      reform_bar2 = extra_partition_list$reform_bar2, betahat = betahat, cov_betahat = cov_betahat,
      dim_coords = object$dim_coords, contrasts = object$contrasts, formula = object$terms,
      xlevels = object$xlevels, diagtol = object$diagtol, family = object$family,
      dispersion_params_val = dispersion_params_val,
      randcov_terms = extra_randcov_list$randcov_terms,
      partition_index_obdata = extra_partition_list$partition_index_obdata
    )

    if (local_list$parallel) {
      # turn of parallel as it is used different in predict
      local_list$parallel <- FALSE
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), loocv_local_glm, loocv_context, se.fit, local_list)
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), loocv_local_glm, loocv_context, se.fit, local_list)
    }
    if (se.fit) {
      cv_predict_val <- vapply(cv_predict_val_list, function(x) x$fit, numeric(1))
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    } else {
      cv_predict_val <- unlist(cv_predict_val_list)
    }
  }

  # loocv computations above happen on the link scale (like the latent w),
  # so map back to the response scale via the family's inverse link before
  # comparing against the observed response y for the fit statistics
  cv_predict_val_invlink <- invlink(cv_predict_val, object$family, object$size)


  cv_predict_error <- y - cv_predict_val_invlink
  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)

  loocv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE
  )

  if (!cv_predict && !se.fit) {
    return(loocv_stats)
  } else {
    loocv_out <- list()
    loocv_out$stats <- loocv_stats


    if (cv_predict) {
      if (type == "link") {
        loocv_out$cv_predict <- cv_predict_val
      } else if (type == "response") {
        loocv_out$cv_predict <- cv_predict_val_invlink
      } else {
        stop("Invalid type argument.", call. = FALSE)
      }
    }

    if (se.fit) {
      loocv_out$se.fit <- as.vector(cv_predict_se)
      if (type == "response" && delta) {
        loocv_out$se.fit <- get_delta_se(cv_predict_val, loocv_out$se.fit, object$family, size)
      }
    }

    return(loocv_out)
  }
}

#' @rdname loocv
#' @method loocv spgautor
#' @order 5
#' @export
loocv.spgautor <- function(object, cv_predict = FALSE, type = c("link", "response"), se.fit = FALSE, delta = FALSE, local, ...) {
  # match type argument so the two display
  type <- match.arg(type)

  if (missing(local)) local <- NULL

  local_list <- get_local_list_prediction(local)

  cov_matrix_val <- covmatrix(object) # already subsets by observed
  X <- model.matrix(object)
  y <- object$y
  cholprods <- get_cholprods_glm(cov_matrix_val, X, y)
  # actually need inverse because of HW blocking
  SigInv <- chol2inv(cholprods$Sig_lowchol)
  SigInv_X <- backsolve(t(cholprods$Sig_lowchol), cholprods$SqrtSigInv_X)

  # find products
  Xt_SigInv_X <- crossprod(X, SigInv_X)
  Xt_SigInv_X_upchol <- base::chol(Xt_SigInv_X) # or Matrix::forceSymmetric()
  cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

  # glm stuff
  dispersion <- as.vector(coef(object, type = "dispersion")) # take class away
  w_linpred <- fitted(object, type = "link") # offset included
  size <- object$size

  # kriging covaries via Sigma (the offset-free latent process) while
  # get_d()/get_D() differentiate the data model; see w_offset_free()
  model_offset <- model.offset(model.frame(object))
  w <- w_offset_free(w_linpred, model_offset)

  # some products
  SigInv_w <- SigInv %*% w
  wX <- cbind(w, X)
  SigInv_wX <- cbind(SigInv_w, SigInv_X)

  # find H stuff
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  Ptheta <- SigInv - SigInv_X %*% wts_beta
  d <- get_d(object$family, w_linpred, y, size, dispersion)
  # and then the gradient vector
  # g <-  d - Ptheta %*% w
  # Next, compute H
  D <- get_D(object$family, w_linpred, y, size, dispersion)
  H <- D - Ptheta
  mHinv <- solve(-H) # chol2inv(chol(Matrix::forceSymmetric(-H))) # solve(-H)

  # parallel stuff
  if (local_list$parallel) {
    cl <- parallel::makeCluster(local_list$ncores)
    cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv_glm,
      Sig = cov_matrix_val,
      SigInv = SigInv, Xmat = X, w = w, wX = wX,
      SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
    )
    cl <- parallel::stopCluster(cl)
  } else {
    cv_predict_val_list <- lapply(seq_len(object$n), get_loocv_glm,
      Sig = cov_matrix_val,
      SigInv = SigInv, Xmat = X, w = as.matrix(w, ncol = 1), wX = wX,
      SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
    )
  }
  cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
  # give each held-out row its own offset back (see the analogous step in
  # loocv.spglm())
  if (!is.null(model_offset)) {
    cv_predict_val <- cv_predict_val + as.vector(model_offset)
  }
  if (se.fit) {
    cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
  }

  cv_predict_val_invlink <- invlink(cv_predict_val, object$family, object$size)

  cv_predict_error <- y - cv_predict_val_invlink
  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)

  loocv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE
  )

  if (!cv_predict && !se.fit) {
    return(loocv_stats)
  } else {
    loocv_out <- list()
    loocv_out$stats <- loocv_stats

    if (cv_predict) {
      if (type == "link") {
        loocv_out$cv_predict <- cv_predict_val
      } else if (type == "response") {
        loocv_out$cv_predict <- cv_predict_val_invlink
      } else {
        stop("Invalid type argument.", call. = FALSE)
      }
    }

    if (se.fit) {
      loocv_out$se.fit <- as.vector(cv_predict_se)
      if (type == "response" && delta) {
        loocv_out$se.fit <- get_delta_se(cv_predict_val, loocv_out$se.fit, object$family, size)
      }
    }

    return(loocv_out)
  }
}

#' Leave-one-out cross validation for a single held-out row (GLM-type big data)
#'
#' @param row The index of the observation to hold out
#' @param context A precomputed context list (built once by \code{loocv.spglm()}
#'   before looping over rows) with the anisotropy-transformed data, model
#'   matrices, fitted covariance/dispersion/random effect parameters, and
#'   \code{randcov_terms}
#' @param se.fit Whether to compute the standard error
#' @param local_list A fully-specified big-data \code{local} list for prediction
#'
#' @return The held-out link-scale fitted value (or a list with \code{fit} and
#'   \code{se.fit} if \code{se.fit} is \code{TRUE})
#'
#' @noRd
loocv_local_glm <- function(row, context, se.fit, local_list) {
  # slice the precomputed context down to the row-removed data instead of
  # calling predict() (which would re-derive all of this from object$formula/
  # object$obdata via model.frame()/model.matrix()/transform_anis() on every
  # single row); level_index_map indexes into group_label, so it must be
  # rebuilt against the row-removed vector rather than reused as-is
  randcov_terms <- lapply(context$randcov_terms, function(term) {
    term$group_label <- term$group_label[-row]
    term$level_index_map <- split(seq_along(term$group_label), term$group_label)
    if (!is.null(term$slope_val)) {
      term$slope_val <- term$slope_val[-row]
    }
    term
  })

  partition_index_obdata <- context$partition_index_obdata
  if (!is.null(partition_index_obdata)) {
    partition_index_obdata$group_label <- partition_index_obdata$group_label[-row]
    partition_index_obdata$level_index_map <- split(
      seq_along(partition_index_obdata$group_label), partition_index_obdata$group_label
    )
  }

  prediction_object <- list(
    se.fit = se.fit, interval = "none", formula = context$formula,
    obdata = context$obdata_pred[-row, , drop = FALSE], xcoord = context$xcoord, ycoord = context$ycoord,
    spcov_params_val = context$spcov_params_val, random = context$random,
    randcov_params_val = context$randcov_params_val, randcov_terms = randcov_terms,
    partition_factor = context$partition_factor, reform_bar2 = context$reform_bar2,
    partition_index_obdata = partition_index_obdata, cov_lowchol = NULL,
    Xmat = context$Xmat_full[-row, , drop = FALSE],
    y = context$y_full[-row],
    betahat = context$betahat, cov_betahat = context$cov_betahat,
    dim_coords = context$dim_coords, contrasts = context$contrasts, local = local_list,
    family = context$family, w = context$w_full[-row],
    model_offset = if (is.null(context$model_offset_full)) NULL else context$model_offset_full[-row],
    size = if (is.null(context$size_full)) NULL else context$size_full[-row],
    dispersion = context$dispersion_params_val, predvar_adjust_ind = TRUE,
    xlevels = context$xlevels, diagtol = context$diagtol, type = "link",
    dist_matrix_full = NULL, partition_vector_full = NULL, cov_vector_full = NULL
  )
  pred <- get_pred_spglm(
    newdata_list = list(row = context$obdata_pred[row, , drop = FALSE], x0 = context$Xmat_full[row, , drop = FALSE]),
    prediction_object = prediction_object
  )

  fit <- pred$fit
  # get_pred_spglm() works with w already offset-adjusted (see the "adjust w"
  # step inside it), so add the held-out row's own offset back in here to
  # return fit on the same (offset-inclusive) link scale as object$y/w
  if (!is.null(context$model_offset_full)) {
    fit <- fit + context$model_offset_full[row]
  }

  if (se.fit) {
    list(fit = fit, se.fit = sqrt(pred$var))
  } else {
    fit
  }
}
