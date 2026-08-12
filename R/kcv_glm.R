#' @param type The scale (\code{response} or \code{link}) of predictions obtained
#' when \code{cv_predict = TRUE} and using \code{spglm()} or \code{spgautor} objects.
#' @param delta A logical indicating whether to return delta method standard errors
#' on the response scale when \code{se.fit = TRUE} and \code{type = "response"}. The default is \code{FALSE}.
#' @rdname kcv
#' @method kcv spglm
#' @order 4
#' @export
kcv.spglm <- function(object, k = 5, cv_predict = FALSE, type = c("link", "response"), se.fit = FALSE, delta = FALSE, local, folds_index, ...) {
  # match type argument so the two display
  type <- match.arg(type)

  if (missing(local)) local <- NULL
  if (missing(folds_index)) folds_index <- NULL
  
  if (is.null(folds_index)) {
    check_kcv_k(k, object$n)

    if (k == object$n) {
      return(loocv(object, cv_predict = cv_predict, type = type, se.fit = se.fit, delta = delta, local = local, ...))
    }
    fold_id <- get_kcv_folds(k, object$n)
  } else {
    check_kcv_folds_index(folds_index, object$n)
    fold_id <- folds_index
  }

  if (is.null(local)) {
    if (object$n > 5000) {
      local <- TRUE
      message("Because the sample size exceeds 5000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun kcv() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }
  local_list <- get_local_list_prediction(local)

  fold_list <- split(seq_len(object$n), fold_id)

  y <- object$y

  cv_predict_val <- numeric(object$n)
  if (se.fit) cv_predict_se <- numeric(object$n)

  if (local_list$method == "all") {
    # exact k-fold CV: build the covariance matrix (on the latent link scale)
    # and its Cholesky-based inverse once, then reuse for every fold via
    # get_kcv_glm()'s closed-form (block Helmert-Wolf-Blocking) update instead
    # of refitting k times
    cov_matrix_val <- covmatrix(object)
    X <- model.matrix(object)
    cholprods <- get_cholprods_glm(cov_matrix_val, X, y)
    SigInv <- chol2inv(cholprods$Sig_lowchol)
    SigInv_X <- backsolve(t(cholprods$Sig_lowchol), cholprods$SqrtSigInv_X)

    # find products
    Xt_SigInv_X <- crossprod(X, SigInv_X)
    Xt_SigInv_X_upchol <- base::chol(Xt_SigInv_X)
    cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

    # glm stuff
    dispersion <- as.vector(coef(object, type = "dispersion"))
    w <- fitted(object, type = "link")
    size <- object$size

    # some products
    SigInv_w <- SigInv %*% w
    wX <- cbind(w, X)
    SigInv_wX <- cbind(SigInv_w, SigInv_X)

    # find H stuff
    wts_beta <- tcrossprod(cov_betahat, SigInv_X)
    Ptheta <- SigInv - SigInv_X %*% wts_beta
    d <- get_d(object$family, w, y, size, dispersion)
    D <- get_D(object$family, w, y, size, dispersion)
    H <- D - Ptheta
    mHinv <- solve(-H)

    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, fold_list, get_kcv_glm,
        Sig = cov_matrix_val,
        SigInv = SigInv, Xmat = X, w = w, wX = wX,
        SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
      )
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(fold_list, get_kcv_glm,
        Sig = cov_matrix_val,
        SigInv = SigInv, Xmat = X, w = as.matrix(w, ncol = 1), wX = wX,
        SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
      )
    }
    for (i in seq_along(fold_list)) {
      fold_rows <- fold_list[[i]]
      cv_predict_val[fold_rows] <- cv_predict_val_list[[i]]$pred
      if (se.fit) cv_predict_se[fold_rows] <- cv_predict_val_list[[i]]$se.fit
    }
  } else {
    # local/big data: refit per fold with the covariance and dispersion
    # parameters held known (see kcv.splm()'s local branch for why this is
    # cheap regardless of n), predicting the now-missing fold on the link
    # scale via predict(), which passes local through to its own
    # nearest-neighbor approximation
    response_name <- all.vars(object$formula)[1]

    spcov_params_val <- coef(object, type = "spcov")
    randcov_params_val <- coef(object, type = "randcov")
    dispersion_params_val <- as.vector(coef(object, type = "dispersion"))

    spcov_initial_val <- do.call(
      spcov_initial,
      c(list(spcov_type = class(spcov_params_val)), as.list(spcov_params_val), list(known = "given"))
    )
    randcov_initial_val <- if (is.null(object$random)) {
      NULL
    } else {
      do.call(randcov_initial, c(as.list(randcov_params_val), list(known = "given")))
    }
    dispersion_initial_val <- dispersion_initial(object$family, dispersion = dispersion_params_val, known = "given")

    for (fold_rows in fold_list) {
      data_train <- object$obdata
      data_train[[response_name]][fold_rows] <- NA

      # xcoord/ycoord must go through do.call() -- see kcv.splm()'s local
      # branch for why a direct call would break
      refit <- do.call("spglm", list(
        formula = object$formula, data = data_train,
        spcov_initial = spcov_initial_val, dispersion_initial = dispersion_initial_val,
        xcoord = object$xcoord, ycoord = object$ycoord, estmethod = object$estmethod,
        anisotropy = object$anisotropy, random = object$random,
        randcov_initial = randcov_initial_val, partition_factor = object$partition_factor,
        local = local
      ))

      # refit$missing_index (ascending) always equals fold_rows (ascending),
      # so predict()'s default no-newdata output lines up positionally
      pred <- predict(refit, type = "link", se.fit = se.fit, local = local, interval = "none")
      if (se.fit) {
        cv_predict_val[fold_rows] <- pred$fit
        cv_predict_se[fold_rows] <- pred$se.fit
      } else {
        cv_predict_val[fold_rows] <- pred
      }
    }
  }

  # kcv computations above happen on the link scale (like the latent w),
  # so map back to the response scale via the family's inverse link before
  # comparing against the observed response y for the fit statistics
  cv_predict_val_invlink <- invlink(cv_predict_val, object$family, object$size)

  cv_predict_error <- y - cv_predict_val_invlink
  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)

  kcv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE
  )

  if (!cv_predict && !se.fit) {
    return(kcv_stats)
  } else {
    kcv_out <- list()
    kcv_out$stats <- kcv_stats

    if (cv_predict) {
      if (type == "link") {
        kcv_out$cv_predict <- cv_predict_val
      } else if (type == "response") {
        kcv_out$cv_predict <- cv_predict_val_invlink
      } else {
        stop("Invalid type argument.", call. = FALSE)
      }
    }

    if (se.fit) {
      kcv_out$se.fit <- as.vector(cv_predict_se)
      if (type == "response" && delta) {
        # object$size (not a possibly-path-local "size") -- the local/refit
        # branch above never defines one
        kcv_out$se.fit <- get_delta_se(cv_predict_val, kcv_out$se.fit, object$family, object$size)
      }
    }

    return(kcv_out)
  }
}

#' @rdname kcv
#' @method kcv spgautor
#' @order 5
#' @export
kcv.spgautor <- function(object, k, cv_predict = FALSE, type = c("link", "response"), se.fit = FALSE, delta = FALSE, local, folds_index, ...) {
  # match type argument so the two display
  type <- match.arg(type)

  if (missing(k)) k <- 5

  if (missing(local)) local <- NULL
  if (missing(folds_index)) folds_index <- NULL

  if (is.null(folds_index)) {
    check_kcv_k(k, object$n)

    if (k == object$n) {
      return(loocv(object, cv_predict = cv_predict, type = type, se.fit = se.fit, delta = delta, local = local, ...))
    }
    fold_id <- get_kcv_folds(k, object$n)
  } else {
    check_kcv_folds_index(folds_index, object$n)
    fold_id <- folds_index
  }

  # spgautor() has no big-data local approximation path -- local is only ever
  # used below for local_list$parallel, matching loocv.spgautor()'s existing
  # convention
  local_list <- get_local_list_prediction(local)

  fold_list <- split(seq_len(object$n), fold_id)

  cov_matrix_val <- covmatrix(object) # already subsets by observed
  X <- model.matrix(object)
  y <- object$y
  cholprods <- get_cholprods_glm(cov_matrix_val, X, y)
  SigInv <- chol2inv(cholprods$Sig_lowchol)
  SigInv_X <- backsolve(t(cholprods$Sig_lowchol), cholprods$SqrtSigInv_X)

  # find products
  Xt_SigInv_X <- crossprod(X, SigInv_X)
  Xt_SigInv_X_upchol <- base::chol(Xt_SigInv_X)
  cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

  # glm stuff
  dispersion <- as.vector(coef(object, type = "dispersion"))
  w <- fitted(object, type = "link")
  size <- object$size

  # some products
  SigInv_w <- SigInv %*% w
  wX <- cbind(w, X)
  SigInv_wX <- cbind(SigInv_w, SigInv_X)

  # find H stuff
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  Ptheta <- SigInv - SigInv_X %*% wts_beta
  d <- get_d(object$family, w, y, size, dispersion)
  D <- get_D(object$family, w, y, size, dispersion)
  H <- D - Ptheta
  mHinv <- solve(-H)

  if (local_list$parallel) {
    cl <- parallel::makeCluster(local_list$ncores)
    cv_predict_val_list <- parallel::parLapply(cl, fold_list, get_kcv_glm,
      Sig = cov_matrix_val,
      SigInv = SigInv, Xmat = X, w = w, wX = wX,
      SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
    )
    cl <- parallel::stopCluster(cl)
  } else {
    cv_predict_val_list <- lapply(fold_list, get_kcv_glm,
      Sig = cov_matrix_val,
      SigInv = SigInv, Xmat = X, w = as.matrix(w, ncol = 1), wX = wX,
      SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
    )
  }

  cv_predict_val <- numeric(object$n)
  if (se.fit) cv_predict_se <- numeric(object$n)
  for (i in seq_along(fold_list)) {
    fold_rows <- fold_list[[i]]
    cv_predict_val[fold_rows] <- cv_predict_val_list[[i]]$pred
    if (se.fit) cv_predict_se[fold_rows] <- cv_predict_val_list[[i]]$se.fit
  }

  cv_predict_val_invlink <- invlink(cv_predict_val, object$family, object$size)

  cv_predict_error <- y - cv_predict_val_invlink
  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)

  kcv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE
  )

  if (!cv_predict && !se.fit) {
    return(kcv_stats)
  } else {
    kcv_out <- list()
    kcv_out$stats <- kcv_stats

    if (cv_predict) {
      if (type == "link") {
        kcv_out$cv_predict <- cv_predict_val
      } else if (type == "response") {
        kcv_out$cv_predict <- cv_predict_val_invlink
      } else {
        stop("Invalid type argument.", call. = FALSE)
      }
    }

    if (se.fit) {
      kcv_out$se.fit <- as.vector(cv_predict_se)
      if (type == "response" && delta) {
        kcv_out$se.fit <- get_delta_se(cv_predict_val, kcv_out$se.fit, object$family, size)
      }
    }

    return(kcv_out)
  }
}
