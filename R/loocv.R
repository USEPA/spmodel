#' Perform leave-one-out cross validation
#'
#' @description Perform leave-one-out cross validation with options for computationally
#'   efficient approximations for big data.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param cv_predict A logical indicating whether the leave-one-out fitted values
#'   should be returned. Defaults to \code{FALSE}. If \code{object} is from [spglm()] or [spgautor()],
#'   the fitted values returned are on the link scale.
#' @param se.fit A logical indicating whether the leave-one-out
#'   prediction standard errors should be returned. Defaults to \code{FALSE}.
#'   If \code{object} is from [spglm()] or [spgautor()],
#'   the standard errors correspond to the fitted values returned on the link scale.
#' @param local A list or logical. If a list, specific list elements described
#'   in [predict.spmodel()] control the big data approximation behavior.
#'   If a logical, \code{TRUE} chooses default list elements for the list version
#'   of \code{local} as specified in [predict.spmodel()]. Defaults to \code{FALSE},
#'   which performs exact computations.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Each observation is held-out from the data set and the remaining data
#'   are used to make a prediction for the held-out observation. This is compared
#'   to the true value of the observation and a mean-squared error is computed
#'   across all observations. The lower the mean squared error, the better the
#'   model fit (according to the leave-one-out criterion).
#'
#' @return If \code{cv_predict = FALSE} and \code{se.fit = FALSE},
#'   a numeric vector indicating the mean-squared-prediction
#'   leave-one-out cross validation error. If \code{cv_predict = TRUE} or \code{se.fit = TRUE},
#'   a list with elements: \code{mspe}, a numeric vector indicating the mean-squared-prediction
#'   leave-one-out cross validation error; \code{cv_predict}, a numeric vector
#'   with leave-one-out predictions for each observation (if \code{cv_predict = TRUE});
#'   and \code{se.fit}, a numeric vector with leave-one-out prediction standard
#'   errors for each observation (if \code{se.fit = TRUE}).
#'
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' loocv(spmod)
#' loocv(spmod, cv_predict = TRUE, se.fit = TRUE)
loocv <- function(object, ...) {
  UseMethod("loocv", object)
}

#' @rdname loocv
#' @method loocv splm
#' @order 2
#' @export
loocv.splm <- function(object, cv_predict = FALSE, se.fit = FALSE, local, ...) {
  if (missing(local)) {
    local <- NULL
  }

  # iid if relevant otherwise pass
  if (inherits(coef(object, type = "spcov"), "none") && is.null(object$random)) {
    return(loocv_iid(object, cv_predict, se.fit, local))
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

  if (local_list$method == "all") {
    # spcov_params_val <- coef(object, type = "spcov")
    # dist_matrix <- spdist(object$obdata, object$xcoord, object$ycoord)
    # randcov_params_val <- coef(object, type = "randcov")
    # if (is.null(object$random)) {
    #   randcov_names <- NULL
    #   randcov_Zs <- NULL
    # } else {
    #   randcov_names <- get_randcov_names(object$random)
    #   randcov_Zs <- get_randcov_Zs(object$obdata, randcov_names)
    # }
    # partition_matrix_val <- partition_matrix(object$partition_factor, object$obdata)
    # cov_matrix_val <- cov_matrix(
    #   spcov_params_val, dist_matrix, randcov_params_val,
    #   randcov_Zs, partition_matrix_val
    # )
    cov_matrix_val <- covmatrix(object)

    # actually need inverse because of HW blocking
    cov_matrixInv_val <- chol2inv(chol(forceSymmetric(cov_matrix_val)))
    model_frame <- model.frame(object)
    X <- model.matrix(object)
    y <- model.response(model_frame)
    yX <- cbind(y, X)
    SigInv_yX <- cov_matrixInv_val %*% yX

    # parallel stuff
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX, se.fit = se.fit
      )
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), get_loocv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX, se.fit = se.fit
      )
    }
    # cv_predict_val <- unlist(cv_predict_val_list)
    cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    if (se.fit) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  } else {
    model_frame <- model.frame(object)
    y <- model.response(model_frame)

    if (local_list$parallel) {
      # turn of parallel as it is used different in predict
      local_list$parallel <- FALSE
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), loocv_local, object, se.fit, local_list)
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), loocv_local, object, se.fit, local_list)
    }
    if (se.fit) {
      cv_predict_val <- vapply(cv_predict_val_list, function(x) x$fit, numeric(1))
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    } else {
      cv_predict_val <- unlist(cv_predict_val_list)
    }
  }
  if (cv_predict) {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), cv_predict = as.vector(cv_predict_val), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), cv_predict = as.vector(cv_predict_val))
    }
  } else {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- mean((cv_predict_val - y)^2)
    }
  }
  cv_output
}

#' @rdname loocv
#' @method loocv spautor
#' @order 3
#' @export
loocv.spautor <- function(object, cv_predict = FALSE, se.fit = FALSE, local, ...) {
  if (missing(local)) {
    local <- NULL
  }

  local_list <- get_local_list_prediction(local)

  # local not used but needed for S3
  # spcov_params_val <- coef(object, type = "spcov")
  # dist_matrix <- object$W
  # randcov_params_val <- coef(object, type = "randcov")
  # if (is.null(object$random)) {
  #   randcov_names <- NULL
  #   randcov_Zs <- NULL
  # } else {
  #   randcov_names <- get_randcov_names(object$random)
  #   randcov_Zs <- get_randcov_Zs(object$data, randcov_names)
  # }
  # partition_matrix_val <- partition_matrix(object$partition_factor, object$data)
  # cov_matrix_val <- cov_matrix(
  #   spcov_params_val, dist_matrix, randcov_params_val,
  #   randcov_Zs, partition_matrix_val, object$M
  # )
  # cov_matrix_obs_val <- cov_matrix_val[object$observed_index, object$observed_index, drop = FALSE]
  cov_matrix_obs_val <- covmatrix(object)

  # actually need inverse because of HW blocking
  cov_matrixInv_obs_val <- chol2inv(chol(forceSymmetric(cov_matrix_obs_val)))
  model_frame <- model.frame(object)
  X <- model.matrix(object)
  y <- model.response(model_frame)
  yX <- cbind(y, X)
  SigInv_yX <- cov_matrixInv_obs_val %*% yX
  # parallel stuff
  if (local_list$parallel) {
    cl <- parallel::makeCluster(local_list$ncores)
    cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv,
      Sig = cov_matrix_obs_val,
      SigInv = cov_matrixInv_obs_val, Xmat = X, y = y, yX = yX,
      SigInv_yX = SigInv_yX, se.fit = se.fit
    )
    cl <- parallel::stopCluster(cl)
  } else {
    cv_predict_val_list <- lapply(seq_len(object$n), get_loocv,
      Sig = cov_matrix_obs_val,
      SigInv = cov_matrixInv_obs_val, Xmat = X, y = y, yX = yX,
      SigInv_yX = SigInv_yX, se.fit = se.fit
    )
  }
  # cv_predict_val <- unlist(cv_predict_val_list)
  cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
  if (se.fit) {
    cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
  }
  if (cv_predict) {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), cv_predict = as.vector(cv_predict_val), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), cv_predict = as.vector(cv_predict_val))
    }
  } else {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- mean((cv_predict_val - y)^2)
    }
  }
  cv_output
}

loocv_local <- function(row, object, se.fit, local_list) {
  newdata <- object$obdata[row, , drop = FALSE]
  object$obdata <- object$obdata[-row, , drop = FALSE]
  predict(object, newdata = newdata, se.fit = se.fit, local = local_list)
}

loocv_iid <- function(object, cv_predict, se.fit, local) {

  # set to FALSE unless it is a list with parallel
  if (is.null(local) || is.logical(local)) local <- FALSE
  local_list <- get_local_list_prediction(local)

  model_frame <- model.frame(object)
  X <- model.matrix(object)
  y <- model.response(model_frame)
  loocv_error <- residuals(object) / (1 - hatvalues(object))
  cv_predict_val <- y - loocv_error

  # parallel stuff
  if (se.fit) {
    total_var <- coef(object, type = "spcov")[["ie"]]
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_se_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv_iid_se,
                                                 vcov(object), Xmat = X, total_var = total_var)
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_se_list <- lapply(seq_len(object$n), get_loocv_iid_se, vcov(object),
                                    Xmat = X, total_var = total_var)
    }
    cv_predict_se <- vapply(cv_predict_se_list, function(x) x$se.fit, numeric(1))
  }

  if (cv_predict) {
    if (se.fit) {
      cv_output <- list(mspe = mean((loocv_error)^2), cv_predict = as.vector(cv_predict_val), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- list(mspe = mean((loocv_error)^2), cv_predict = as.vector(cv_predict_val))
    }
  } else {
    if (se.fit) {
      cv_output <- list(mspe = mean((loocv_error)^2), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- mean((loocv_error)^2)
    }
  }
  cv_output
}
