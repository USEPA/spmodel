#' Perform leave-one-out cross validation
#'
#' @description Perform leave-one-out cross validation with options for computationally
#'   efficient approximations for big data.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param cv_predict A logical indicating whether the leave-one-out fitted values
#'   should be returned. Defaults to \code{FALSE}.
#' @param local A list or logical. If a list, specific list elements described
#'   in [predict.spmod()] control the big data approximation behavior.
#'   If a logical, \code{TRUE} chooses default list elements for the list version
#'   of \code{local} as specified in [predict.spmod()]. Defaults to \code{FALSE},
#'   which performs exact computations.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Each observation is held-out from the data set and the remaining data
#'   are used to make a prediction for the held-out observation. This is compared
#'   to the true value of the observation and a mean-squared error is computed
#'   across all observations. The lower the mean squared error, the better the
#'   model fit (according to the leave-one-out criterion).
#'
#' @return If \code{cv_predict = FALSE}, a numeric vector indicating the mean-squared-prediction
#'   leave-one-out cross validation error. If \code{cv_predict = TRUE}, a list with
#'   two elements: \code{mspe}, a numeric vector indicating the mean-squared-prediction
#'   leave-one-out cross validation error; and \code{cv_predict}, a numeric vector
#'   with leave-one-out predictions for each observation.
#'
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' loocv(spmod)
#' loocv(spmod, cv_predict = TRUE)
loocv <- function(object, ...) {
  UseMethod("loocv", object)
}

#' @rdname loocv
#' @method loocv spmod
#' @export
loocv.spmod <- function(object, cv_predict = FALSE, local, ...) {
  if (missing(local)) {
    local <- NULL
  }

  switch(object$fn,
    "splm" = loocv_splm(object, cv_predict, local, ...),
    "spautor" = loocv_spautor(object, cv_predict, local, ...)
  )
}

loocv_splm <- function(object, cv_predict = FALSE, local, ...) {

  # local prediction list
  # browser()

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
    spcov_params_val <- coef(object, type = "spcov")
    dist_matrix <- spdist(object$obdata, object$xcoord, object$ycoord)
    randcov_params_val <- coef(object, type = "randcov")
    if (is.null(object$random)) {
      randcov_names <- NULL
      randcov_Zs <- NULL
    } else {
      randcov_names <- get_randcov_names(object$random)
      randcov_Zs <- get_randcov_Zs(object$obdata, randcov_names)
    }
    partition_matrix_val <- partition_matrix(object$partition_factor, object$obdata)
    cov_matrix_val <- cov_matrix(
      spcov_params_val, dist_matrix, randcov_params_val,
      randcov_Zs, partition_matrix_val
    )

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
        SigInv_yX = SigInv_yX
      )
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), get_loocv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX
      )
    }
    cv_predict_val <- unlist(cv_predict_val_list)
  } else {
    model_frame <- model.frame(object)
    y <- model.response(model_frame)

    if (local_list$parallel) {
      # turn of parallel as it is used different in predict
      local_list$parallel <- FALSE
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), loocv_local, object, local_list)
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), loocv_local, object, local_list)
    }
    cv_predict_val <- unlist(cv_predict_val_list)
    # cv_predict_val <- do.call("rbind", cv_predict_val)
  }
  if (cv_predict) {
    cv_output <- list(mspe = mean((cv_predict_val - y)^2), cv_predict = as.vector(cv_predict_val))
  } else {
    cv_output <- mean((cv_predict_val - y)^2)
  }
  cv_output
}

loocv_spautor <- function(object, cv_predict = FALSE, local, ...) {
  local_list <- get_local_list_prediction(local)

  # local not used but needed for S3
  spcov_params_val <- coef(object, type = "spcov")
  dist_matrix <- object$W
  randcov_params_val <- coef(object, type = "randcov")
  if (is.null(object$random)) {
    randcov_names <- NULL
    randcov_Zs <- NULL
  } else {
    randcov_names <- get_randcov_names(object$random)
    randcov_Zs <- get_randcov_Zs(object$data, randcov_names)
  }
  partition_matrix_val <- partition_matrix(object$partition_factor, object$data)
  cov_matrix_val <- cov_matrix(
    spcov_params_val, dist_matrix, randcov_params_val,
    randcov_Zs, partition_matrix_val, object$M
  )
  cov_matrix_obs_val <- cov_matrix_val[object$observed_index, object$observed_index, drop = FALSE]

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
      Sig = cov_matrix_val,
      SigInv = cov_matrixInv_obs_val, Xmat = X, y = y, yX = yX,
      SigInv_yX = SigInv_yX
    )
    cl <- parallel::stopCluster(cl)
  } else {
    cv_predict_val_list <- lapply(seq_len(object$n), get_loocv,
      Sig = cov_matrix_val,
      SigInv = cov_matrixInv_obs_val, Xmat = X, y = y, yX = yX,
      SigInv_yX = SigInv_yX
    )
  }
  cv_predict_val <- unlist(cv_predict_val_list)
  if (cv_predict) {
    cv_output <- list(mspe = mean((cv_predict_val - y)^2), cv_predict = as.vector(cv_predict_val))
  } else {
    cv_output <- mean((cv_predict_val - y)^2)
  }
  cv_output
}

#   # maybe put cv fitted in fitted?
loocv_local <- function(row, object, local_list) {
  newdata <- object$obdata[row, , drop = FALSE]
  object$obdata <- object$obdata[-row, , drop = FALSE]
  predict(object, newdata = newdata, local = local_list)
}
