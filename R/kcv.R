#' Perform k-fold cross validation
#'
#' @description Perform k-fold cross validation with options for computationally
#'   efficient approximations for big data. Generalizes [loocv()] (leave-one-out
#'   cross validation) to leaving out \code{k} folds of (approximately) equal size
#'   instead of single observations.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param k The number of folds. Must be a whole number at least 2 and no more
#'   than the sample size (the number of non-missing observations in \code{object}).
#'   The default is \code{5}. If \code{k} equals the sample size, \code{kcv()}
#'   is equivalent to (and calls) [loocv()].
#' @param cv_predict A logical indicating whether the k-fold cross validation fitted values
#'   should be returned. Defaults to \code{FALSE}. If \code{object} is from [spglm()] or [spgautor()],
#'   the fitted values returned are on the link scale.
#' @param se.fit A logical indicating whether the k-fold cross validation
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
#' @details Observations are randomly partitioned into \code{k} folds of
#'   (approximately) equal size. Each fold is held out from the data set in turn
#'   and the remaining data are used to make predictions for the held-out fold. This is
#'   compared to the true values of the held-out observations and several fit
#'   statistics are (sometimes optionally) computed: bias, mean-squared-prediction error (MSPE),
#'   root-mean-squared-prediction error (RMSPE), and the squared correlation
#'   (cor2) between the observed data and k-fold cross validation predictions
#'   (regarded as a prediction version of r-squared appropriate for comparing
#'   across spatial and nonspatial models), , and 
#'   prediction interval coverage (cover.XX). Generally,
#'   bias should be near zero and prediction interval coverage at the
#'   intended level for well-fitting models. The lower the MSPE and RMSPE, the better the model
#'   fit (according to the k-fold cross validation criterion). The higher the
#'   cor2, the better the model fit (according to the k-fold cross validation
#'   criterion). cor2 and cover.XX are not returned when \code{object} was fit using
#'   \code{spglm()} or \code{spgautor()} because we do not observe the underlying latent mean.
#'
#'   When \code{object} is from \code{splm()} or \code{spautor()}, setting
#'   \code{interval = "prediction"} additionally reports the empirical coverage
#'   of the \code{level} (e.g. 95\%) k-fold cross validation prediction interval --
#'   the proportion of held-out observations whose true value falls within
#'   \code{fit +/- qnorm(1 - (1 - level) / 2) * se.fit} (the same normal-quantile
#'   interval [predict.spmodel()] uses by default). This is only available for
#'   \code{splm()}/\code{spautor()} objects, since \code{spglm()}/\code{spgautor()}
#'   have no observed-scale latent mean to compare against.
#'
#' @return If \code{cv_predict = FALSE} and \code{se.fit = FALSE},
#'   a fit statistics tibble (with bias, MSPE, RMSPE, and cor2; see Details).
#'   If \code{cv_predict = TRUE} or \code{se.fit = TRUE},
#'   a list with elements: \code{stats}, a fit statistics tibble
#'   (with bias, MSPE, RMSPE, and cor2; see Details); \code{cv_predict}, a numeric vector
#'   with k-fold cross validation predictions for each observation (if \code{cv_predict = TRUE});
#'   and \code{se.fit}, a numeric vector with k-fold cross validation prediction standard
#'   errors for each observation (if \code{se.fit = TRUE}). When \code{object} is from
#'   \code{splm()} or \code{spautor()} and \code{interval = "prediction"}, the fit
#'   statistics tibble also has a \code{cover.XX} column (e.g. \code{cover.95}
#'   for \code{level = 0.95}; see Details).
#'
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' kcv(spmod)
#' kcv(spmod, k = 10, cv_predict = TRUE, se.fit = TRUE)
kcv <- function(object, ...) {
  UseMethod("kcv", object)
}

#' @param interval Whether to also report empirical k-fold cross validation
#'   prediction interval coverage in the returned fit statistics. \code{"none"}
#'   (the default) omits it; \code{"prediction"} reports it (see Details). Only
#'   available for \code{splm()}/\code{spautor()} objects.
#' @param level The prediction interval level (e.g. 0.95) used to compute
#'   prediction interval coverage when \code{interval = "prediction"}. Ignored otherwise. The
#'   default is \code{0.95}.
#' @rdname kcv
#' @method kcv splm
#' @order 2
#' @export
kcv.splm <- function(object, k = 5, cv_predict = FALSE, se.fit = FALSE, local, interval = c("none", "prediction"), level = 0.95, ...) {
  interval <- match.arg(interval)

  check_kcv_k(k, object$n)

  if (missing(local)) local <- NULL

  # k = n is exactly leave-one-out -- delegate rather than duplicating loocv()'s
  # own fast paths (e.g. the closed-form iid shortcut)
  if (k == object$n) {
    return(loocv(object, cv_predict = cv_predict, se.fit = se.fit, local = local, interval = interval, level = level, ...))
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

  fold_id <- get_kcv_folds(k, object$n)
  fold_list <- split(seq_len(object$n), fold_id)

  model_frame <- model.frame(object)
  y <- model.response(model_frame)

  # interval = "prediction" needs se.fit internally to build the interval even
  # when the caller didn't request se.fit in the return value
  se_needed <- se.fit || interval == "prediction"

  cv_predict_val <- numeric(object$n)
  if (se_needed) cv_predict_se <- numeric(object$n)

  if (local_list$method == "all") {
    # exact k-fold CV: build the full n x n covariance matrix and its inverse
    # once, then reuse them for every fold via get_kcv()'s closed-form
    # leave-fold-out (block Helmert-Wolf-Blocking) update instead of refitting
    # the model k times
    cov_matrix_val <- covmatrix(object)
    cov_matrixInv_val <- chol2inv(chol(forceSymmetric(cov_matrix_val)))
    X <- model.matrix(object)
    yX <- cbind(y, X)
    SigInv_yX <- cov_matrixInv_val %*% yX

    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, fold_list, get_kcv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX, se.fit = se_needed
      )
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(fold_list, get_kcv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX, se.fit = se_needed
      )
    }
    for (i in seq_along(fold_list)) {
      fold_rows <- fold_list[[i]]
      cv_predict_val[fold_rows] <- cv_predict_val_list[[i]]$pred
      if (se_needed) cv_predict_se[fold_rows] <- cv_predict_val_list[[i]]$se.fit
    }
  } else {
    # local/big data: avoid ever forming the full n x n covariance matrix.
    # Instead, for each fold, refit with the covariance parameters held known
    # (spcov_initial()/randcov_initial() with known = "given" route fitting to
    # the non-iterative use_gloglik_known() path -- no optim() call, so this is
    # cheap regardless of n) and predict the now-missing fold via predict(),
    # which passes local through to its own nearest-neighbor approximation.
    # betahat is genuinely re-estimated per fold either way (unlike
    # loocv_local(), which keeps betahat fixed -- a negligible approximation
    # for one held-out row, but not for a fold of this size).
    response_name <- all.vars(object$formula)[1]

    spcov_params_val <- coef(object, type = "spcov")
    randcov_params_val <- coef(object, type = "randcov")
    spcov_initial_val <- do.call(
      spcov_initial,
      c(list(spcov_type = class(spcov_params_val)), as.list(spcov_params_val), list(known = "given"))
    )
    randcov_initial_val <- if (is.null(object$random)) {
      NULL
    } else {
      do.call(randcov_initial, c(as.list(randcov_params_val), list(known = "given")))
    }

    for (fold_rows in fold_list) {
      data_train <- object$obdata
      data_train[[response_name]][fold_rows] <- NA

      # splm()'s xcoord/ycoord use non-standard evaluation (substitute()) to
      # support bare unquoted column names, so object$xcoord/object$ycoord
      # (already plain strings) must be passed via do.call() -- a direct call
      # would capture the expression "object$xcoord" instead of its value
      refit <- do.call("splm", list(
        formula = object$formula, data = data_train, spcov_initial = spcov_initial_val,
        xcoord = object$xcoord, ycoord = object$ycoord, estmethod = object$estmethod,
        anisotropy = object$anisotropy, random = object$random,
        randcov_initial = randcov_initial_val, partition_factor = object$partition_factor,
        local = local, ddf = "asymptotic"
      ))

      # refit$missing_index (ascending, from which(is.na(...))) always equals
      # fold_rows (ascending, from split()), so predict()'s default no-newdata
      # output lines up positionally with fold_rows with no reordering needed
      pred <- predict(refit, se.fit = se_needed, local = local, interval = "none")
      if (se_needed) {
        cv_predict_val[fold_rows] <- pred$fit
        cv_predict_se[fold_rows] <- pred$se.fit
      } else {
        cv_predict_val[fold_rows] <- pred
      }
    }
  }

  cv_predict_error <- y - cv_predict_val
  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)
  cor2 <- cor(cv_predict_val, y)^2

  kcv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE,
    cor2 = cor2
  )

  if (interval == "prediction") {
    kcv_stats <- add_cv_coverage(kcv_stats, y, cv_predict_val, cv_predict_se, level)
  }

  if (!cv_predict && !se.fit) {
    return(kcv_stats)
  } else {
    kcv_out <- list()
    kcv_out$stats <- kcv_stats

    if (cv_predict) {
      kcv_out$cv_predict <- cv_predict_val
    }

    if (se.fit) {
      kcv_out$se.fit <- as.vector(cv_predict_se)
    }
    return(kcv_out)
  }
}

#' @rdname kcv
#' @method kcv spautor
#' @order 3
#' @export
kcv.spautor <- function(object, k = 5, cv_predict = FALSE, se.fit = FALSE, local, interval = c("none", "prediction"), level = 0.95, ...) {
  interval <- match.arg(interval)

  check_kcv_k(k, object$n)

  if (missing(local)) local <- NULL

  if (k == object$n) {
    return(loocv(object, cv_predict = cv_predict, se.fit = se.fit, local = local, interval = interval, level = level, ...))
  }

  # spautor() has no big-data local approximation path (areal neighborhood
  # structures are typically much smaller than point-referenced data sets) --
  # local is only ever used below for local_list$parallel, matching
  # loocv.spautor()'s existing convention
  local_list <- get_local_list_prediction(local)
  # interval = "prediction" needs se.fit internally to build the interval even
  # when the caller didn't request se.fit in the return value
  se_needed <- se.fit || interval == "prediction"

  fold_id <- get_kcv_folds(k, object$n)
  fold_list <- split(seq_len(object$n), fold_id)

  cov_matrix_obs_val <- covmatrix(object)
  cov_matrixInv_obs_val <- chol2inv(chol(forceSymmetric(cov_matrix_obs_val)))
  model_frame <- model.frame(object)
  X <- model.matrix(object)
  y <- model.response(model_frame)
  yX <- cbind(y, X)
  SigInv_yX <- cov_matrixInv_obs_val %*% yX

  if (local_list$parallel) {
    cl <- parallel::makeCluster(local_list$ncores)
    cv_predict_val_list <- parallel::parLapply(cl, fold_list, get_kcv,
      Sig = cov_matrix_obs_val,
      SigInv = cov_matrixInv_obs_val, Xmat = X, y = y, yX = yX,
      SigInv_yX = SigInv_yX, se.fit = se_needed
    )
    cl <- parallel::stopCluster(cl)
  } else {
    cv_predict_val_list <- lapply(fold_list, get_kcv,
      Sig = cov_matrix_obs_val,
      SigInv = cov_matrixInv_obs_val, Xmat = X, y = y, yX = yX,
      SigInv_yX = SigInv_yX, se.fit = se_needed
    )
  }

  cv_predict_val <- numeric(object$n)
  if (se_needed) cv_predict_se <- numeric(object$n)
  for (i in seq_along(fold_list)) {
    fold_rows <- fold_list[[i]]
    cv_predict_val[fold_rows] <- cv_predict_val_list[[i]]$pred
    if (se_needed) cv_predict_se[fold_rows] <- cv_predict_val_list[[i]]$se.fit
  }

  cv_predict_error <- y - cv_predict_val
  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)
  cor2 <- cor(cv_predict_val, y)^2

  kcv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE,
    cor2 = cor2
  )

  if (interval == "prediction") {
    kcv_stats <- add_cv_coverage(kcv_stats, y, cv_predict_val, cv_predict_se, level)
  }

  if (!cv_predict && !se.fit) {
    return(kcv_stats)
  } else {
    kcv_out <- list()
    kcv_out$stats <- kcv_stats

    if (cv_predict) {
      kcv_out$cv_predict <- cv_predict_val
    }

    if (se.fit) {
      kcv_out$se.fit <- as.vector(cv_predict_se)
    }
    return(kcv_out)
  }
}

#' Validate the \code{k} argument to \code{kcv()}
#'
#' @param k The user-supplied \code{k} argument
#' @param n The fitted model's sample size
#'
#' @return Invisibly \code{NULL}; called for its error-checking side effect
#'
#' @noRd
check_kcv_k <- function(k, n) {
  if (!is.numeric(k) || length(k) != 1 || k != round(k) || k < 2) {
    stop("k must be a single whole number greater than or equal to 2.", call. = FALSE)
  }
  if (k > n) {
    stop("k cannot exceed the number of observations.", call. = FALSE)
  }
  invisible(NULL)
}

#' Randomly assign observations to (approximately) equally-sized folds
#'
#' Same construction as \code{\link{get_training_list}()}'s \code{"cv"} branch
#' (\code{R/decorrelate.R}), reused here so \code{kcv()}'s fold sizes differ by
#' at most 1.
#'
#' @param k The number of folds
#' @param n The sample size
#'
#' @return A length-\code{n} vector assigning each observation to one of
#'   \code{seq_len(k)}
#'
#' @noRd
get_kcv_folds <- function(k, n) {
  sample(rep(seq_len(k), length.out = n))
}
