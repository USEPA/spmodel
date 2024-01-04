#' @rdname loocv
#' @method loocv spglm
#' @order 4
#' @export
loocv.spglm <- function(object, cv_predict = FALSE, se.fit = FALSE, local, ...) {
  if (missing(local)) {
    local <- NULL
  }

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
    # and then the gradient vector
    # g <-  d - Ptheta %*% w
    # Next, compute H
    D <- get_D(object$family, w, y, size, dispersion)
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
    if (se.fit) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  } else {

    # get w for later
    w <- fitted(object, type = "link")

    extra_randcov_list <- get_extra_randcov_list(object, object$obdata, newdata = object$obdata)
    extra_partition_list <- get_extra_partition_list(object, object$obdata, newdata = object$obdata)

    if (local_list$parallel) {
      # turn of parallel as it is used different in predict
      local_list$parallel <- FALSE
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), loocv_local_glm, object, se.fit, local_list, extra_randcov_list = extra_randcov_list, extra_partition_list = extra_partition_list)
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), loocv_local_glm, object, se.fit, local_list,
                                    extra_randcov_list = extra_randcov_list, extra_partition_list = extra_partition_list)
    }
    if (se.fit) {
      cv_predict_val <- vapply(cv_predict_val_list, function(x) x$fit, numeric(1))
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    } else {
      cv_predict_val <- unlist(cv_predict_val_list)
    }
  }

  cv_predict_val_invlink <- invlink(cv_predict_val, object$family, object$size)

  if (cv_predict) {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), cv_predict = as.vector(cv_predict_val), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), cv_predict = as.vector(cv_predict_val))
    }
  } else {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- mean((cv_predict_val_invlink - y)^2)
    }
  }
  cv_output
}

#' @rdname loocv
#' @method loocv spgautor
#' @order 5
#' @export
loocv.spgautor <- function(object, cv_predict = FALSE, se.fit = FALSE, local, ...) {
  if (missing(local)) {
    local <- NULL
  }

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
  # and then the gradient vector
  # g <-  d - Ptheta %*% w
  # Next, compute H
  D <- get_D(object$family, w, y, size, dispersion)
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
  if (se.fit) {
    cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
  }

  cv_predict_val_invlink <- invlink(cv_predict_val, object$family, object$size)

  if (cv_predict) {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), cv_predict = as.vector(cv_predict_val), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), cv_predict = as.vector(cv_predict_val))
    }
  } else {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- mean((cv_predict_val_invlink - y)^2)
    }
  }
  cv_output
}

loocv_local_glm <- function(row, object, se.fit, local_list,
                            extra_randcov_list = NULL, extra_partition_list = NULL) {
  object$fitted$link <- object$fitted$link[-row] # w needs to be subset
  newdata <- object$obdata[row, , drop = FALSE]
  object$obdata <- object$obdata[-row, , drop = FALSE]
  # this is all so the randcov and partition steps are not repeated for each iteration
  if (!is.null(extra_randcov_list)) {
    extra_randcov_list$Z_index_obdata_list <- lapply(extra_randcov_list$Z_index_obdata_list,
                                                     function(x) {
                                                       x$reform_bar2_vals <- x$reform_bar2_vals[-row]
                                                       x
                                                     })
  }

  if (!is.null(extra_partition_list)) {
    extra_partition_list$partition_index_obdata$reform_bar2_vals <- extra_partition_list$partition_index_obdata$reform_bar2_vals[-row]
  }

  predict(object, newdata = newdata, se.fit = se.fit, local = local_list,
          extra_randcov_list = extra_randcov_list,
          extra_partition_list = extra_partition_list
  )

}
