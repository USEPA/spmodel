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
#'   to the true value of the observation and several fit statistics are (sometimes optionally) computed:
#'   bias, mean-squared-prediction error (MSPE), root-mean-squared-prediction
#'   error (RMSPE), and the squared correlation (cor2) between the observed data
#'   and leave-one-out predictions (regarded as a prediction version of r-squared
#'   appropriate for comparing across spatial and nonspatial models), and 
#'   prediction interval coverage (cover.XX). Generally,
#'   bias should be near zero and prediction interval coverage at the
#'   intended level for well-fitting models. The lower the MSPE and RMSPE,
#'   the better the model fit (according to the leave-out-out criterion).
#'   The higher the cor2, the better the model fit (according to the leave-out-out
#'   criterion). cor2 and cover.XX are not returned when \code{object} was fit using
#'   \code{spglm()} or \code{spgautor()} because we do not observe the underlying latent mean.
#'
#' @return If \code{cv_predict = FALSE} and \code{se.fit = FALSE},
#'   a fit statistics tibble (with bias, MSPE, RMSPE, and cor2; see Details).
#'   If \code{cv_predict = TRUE} or \code{se.fit = TRUE},
#'   a list with elements: \code{stats}, a fit statistics tibble
#'   (with bias, MSPE, RMSPE, and cor2; see Details); \code{cv_predict}, a numeric vector
#'   with leave-one-out predictions for each observation (if \code{cv_predict = TRUE});
#'   and \code{se.fit}, a numeric vector with leave-one-out prediction standard
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
#' loocv(spmod)
#' loocv(spmod, cv_predict = TRUE, se.fit = TRUE)
loocv <- function(object, ...) {
  UseMethod("loocv", object)
}

#' @param interval Whether to also report empirical leave-one-out prediction
#'   interval coverage in the returned fit statistics. \code{"none"} (the
#'   default) omits it; \code{"prediction"} reports it (see Details). Only
#'   available for \code{splm()}/\code{spautor()} objects.
#' @param level The prediction interval level (e.g. 0.95) used to compute
#'   prediction interval coverage when \code{interval = "prediction"}. Ignored otherwise. The
#'   default is \code{0.95}.
#' @rdname loocv
#' @method loocv splm
#' @order 2
#' @export
loocv.splm <- function(object, cv_predict = FALSE, se.fit = FALSE, local, interval = c("none", "prediction"), level = 0.95, ...) {
  interval <- match.arg(interval)
  if (missing(local)) local <- NULL

  # when there is no spatial dependence or random effects, leave-one-out
  # predictions reduce to a closed-form formula (PRESS residuals; see
  # loocv_iid()) that avoids inverting an n x n covariance matrix or looping
  # over rows, so take that fast path whenever it applies
  # iid if relevant otherwise pass
  if (inherits(coef(object, type = "spcov"), c("none", "ie")) && is.null(object$random)) {
    return(loocv_iid(object, cv_predict, se.fit, local, interval, level))
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
  # interval = "prediction" needs se.fit internally to build the interval even
  # when the caller didn't request se.fit in the return value
  se_needed <- se.fit || interval == "prediction"

  if (local_list$method == "all") {
    # exact LOO: build the full n x n covariance matrix and its inverse once,
    # then reuse them for every held-out row via get_loocv()'s closed-form
    # leave-one-out update instead of refitting the model n times
    cov_matrix_val <- covmatrix(object)

    # actually need inverse because of HW blocking
    cov_matrixInv_val <- chol2inv(chol(forceSymmetric(cov_matrix_val)))
    model_frame <- model.frame(object)
    X <- model.matrix(object)
    y <- model.response(model_frame)
    yX <- cbind(y, X)
    # binding y and X together lets a single matrix product below produce both
    # SigInv %*% y and SigInv %*% X, which get_loocv() needs per row
    SigInv_yX <- cov_matrixInv_val %*% yX

    # parallel stuff
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX, se.fit = se_needed
      )
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), get_loocv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX, se.fit = se_needed
      )
    }
    cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    if (se_needed) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  } else {
    # precompute everything that is invariant across held-out rows -- the
    # anisotropy-transformed coordinates, spatial/random covariance params,
    # design matrix, response, offset, betahat, and cov_betahat -- so the
    # per-row worker (loocv_local()) can slice arrays directly instead of
    # calling predict(), which previously re-derived every one of these from
    # object$formula/object$obdata (via model.frame()/model.matrix()/
    # transform_anis()/replace_newdata()) from scratch on every single row.
    # None of these actually depend on which row is held out (betahat and
    # cov_betahat intentionally stay fixed at their full-data estimates for
    # the local approximation, matching the exact -- local = FALSE -- path's
    # existing convention of not refitting spcov/randcov params per row).
    spcov_params_val <- coef(object, type = "spcov")
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

    model_frame <- model.frame(object)
    y <- model.response(model_frame)
    Xmat_full <- model.matrix(object)
    offset_full <- model.offset(model_frame)
    betahat <- coefficients(object)
    cov_betahat <- vcov(object)

    extra_randcov_list <- get_extra_randcov_list(object, object$obdata, newdata = object$obdata)
    extra_partition_list <- get_extra_partition_list(object, object$obdata, newdata = object$obdata)


    # add dummy level if necessary (holding out a row can make it the only
    # observation with a particular factor level, so the level-removed data
    # passed to loocv_local() may need a placeholder level available so
    # model.matrix() doesn't error when that observation is later "predicted")
    ## random effect
    if (!is.null(extra_randcov_list$randcov_terms)) {
      for (x in names(extra_randcov_list$randcov_terms)) {
        val <- extra_randcov_list$randcov_terms[[x]]$xlev[[1]]
        if (!"...this_is_a_new_level..." %in% val) {
          extra_randcov_list$randcov_terms[[x]]$xlev[[1]] <- c(val, "...this_is_a_new_level...")
        }
      }
    }
    ## partition factor
    if (!is.null(extra_partition_list$partition_index_obdata)) {
      val <- extra_partition_list$partition_index_obdata$xlev[[1]]
      if (!"...this_is_a_new_level..." %in% val) {
        extra_partition_list$partition_index_obdata$xlev[[1]] <- c(val, "...this_is_a_new_level...")
      }
    }

    loocv_context <- list(
      obdata_pred = obdata_pred, Xmat_full = Xmat_full, y_full = y, offset_full = offset_full,
      xcoord = xcoord, ycoord = ycoord, spcov_params_val = spcov_params_val, random = object$random,
      randcov_params_val = randcov_params_val, partition_factor = object$partition_factor,
      reform_bar2 = extra_partition_list$reform_bar2, betahat = betahat, cov_betahat = cov_betahat,
      dim_coords = object$dim_coords, contrasts = object$contrasts, formula = object$terms,
      xlevels = object$xlevels, diagtol = object$diagtol,
      randcov_terms = extra_randcov_list$randcov_terms,
      partition_index_obdata = extra_partition_list$partition_index_obdata
    )

    if (local_list$parallel) {
      # turn of parallel as it is used different in predict
      local_list$parallel <- FALSE
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), loocv_local, loocv_context, se_needed, local_list)
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), loocv_local, loocv_context, se_needed, local_list)
    }
    if (se_needed) {
      cv_predict_val <- vapply(cv_predict_val_list, function(x) x$fit, numeric(1))
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    } else {
      cv_predict_val <- unlist(cv_predict_val_list)
    }
  }

  cv_predict_error <- y - cv_predict_val
  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)
  cor2 <- cor(cv_predict_val, y)^2

  loocv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE,
    cor2 = cor2
  )

  if (interval == "prediction") {
    loocv_stats <- add_cv_coverage(loocv_stats, y, cv_predict_val, cv_predict_se, level)
  }

  if (!cv_predict && !se.fit) {
    return(loocv_stats)
  } else {
    loocv_out <- list()
    loocv_out$stats <- loocv_stats

    if (cv_predict) {
      loocv_out$cv_predict <- cv_predict_val
    }

    if (se.fit) {
      loocv_out$se.fit <- as.vector(cv_predict_se)
    }
    return(loocv_out)
  }
}

#' @rdname loocv
#' @method loocv spautor
#' @order 3
#' @export
loocv.spautor <- function(object, cv_predict = FALSE, se.fit = FALSE, local, interval = c("none", "prediction"), level = 0.95, ...) {
  interval <- match.arg(interval)
  if (missing(local)) local <- NULL

  local_list <- get_local_list_prediction(local)
  # interval = "prediction" needs se.fit internally to build the interval even
  # when the caller didn't request se.fit in the return value
  se_needed <- se.fit || interval == "prediction"

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
      SigInv_yX = SigInv_yX, se.fit = se_needed
    )
    cl <- parallel::stopCluster(cl)
  } else {
    cv_predict_val_list <- lapply(seq_len(object$n), get_loocv,
      Sig = cov_matrix_obs_val,
      SigInv = cov_matrixInv_obs_val, Xmat = X, y = y, yX = yX,
      SigInv_yX = SigInv_yX, se.fit = se_needed
    )
  }
  cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
  if (se_needed) {
    cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
  }

  cv_predict_error <- y - cv_predict_val
  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)
  cor2 <- cor(cv_predict_val, y)^2

  loocv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE,
    cor2 = cor2
  )

  if (interval == "prediction") {
    loocv_stats <- add_cv_coverage(loocv_stats, y, cv_predict_val, cv_predict_se, level)
  }

  if (!cv_predict && !se.fit) {
    return(loocv_stats)
  } else {
    loocv_out <- list()
    loocv_out$stats <- loocv_stats

    if (cv_predict) {
      loocv_out$cv_predict <- cv_predict_val
    }

    if (se.fit) {
      loocv_out$se.fit <- as.vector(cv_predict_se)
    }
    return(loocv_out)
  }
}

#' Leave-one-out cross validation for a single held-out row (\code{splm()} big data)
#'
#' @param row The index of the observation to hold out
#' @param context A precomputed context list (built once by \code{loocv.splm()}
#'   before looping over rows) with the anisotropy-transformed data, model
#'   matrices, fitted covariance/random effect parameters, and \code{randcov_terms}
#' @param se.fit Whether to compute the standard error
#' @param local_list A fully-specified big-data \code{local} list for prediction
#'
#' @return The held-out fitted value (or a list with \code{fit} and \code{se.fit}
#'   if \code{se.fit} is \code{TRUE})
#'
#' @noRd
loocv_local <- function(row, context, se.fit, local_list) {
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
    offset = if (is.null(context$offset_full)) NULL else context$offset_full[-row],
    betahat = context$betahat, cov_betahat = context$cov_betahat,
    dim_coords = context$dim_coords, contrasts = context$contrasts, local = local_list,
    xlevels = context$xlevels, diagtol = context$diagtol, type = "response",
    dist_matrix_full = NULL, partition_vector_full = NULL, cov_vector_full = NULL
  )
  pred <- get_pred_splm(
    newdata_list = list(row = context$obdata_pred[row, , drop = FALSE], x0 = context$Xmat_full[row, , drop = FALSE]),
    prediction_object = prediction_object
  )

  fit <- pred$fit
  # get_pred_splm() computed fit on the offset-removed scale (offset was
  # subtracted out of y before forming residuals), so add the held-out row's
  # own offset back in here, mirroring predict.splm()'s "apply offset" step
  if (!is.null(context$offset_full)) {
    fit <- fit + context$offset_full[row]
  }

  if (se.fit) {
    list(fit = fit, se.fit = sqrt(pred$var))
  } else {
    fit
  }
}

#' Leave-one-out cross validation for an \code{splm()} model with iid errors
#'
#' @param object A fitted model object from [splm()]
#' @param cv_predict Whether to return the leave-one-out fitted values
#' @param se.fit Whether to return the leave-one-out prediction standard errors
#' @param local A list or logical controlling the big data approximation for
#'   parallelizing the (independent, so embarrassingly parallel) standard error calculation
#' @param interval See \code{loocv.splm()}'s \code{interval} argument
#' @param level See \code{loocv.splm()}'s \code{level} argument
#'
#' @return The same value as \code{loocv.splm()}, computed via the classical
#'   leave-one-out identity \eqn{y_i - \hat{y}_i^{(-i)} = e_i / (1 - h_{ii})},
#'   since with no spatial dependence or random effects each observation's
#'   leave-one-out residual is a simple function of its ordinary residual and
#'   leverage rather than requiring a full covariance-matrix update
#'
#' @noRd
loocv_iid <- function(object, cv_predict, se.fit, local, interval = "none", level = 0.95) {
  # set to FALSE unless it is a list with parallel
  if (is.null(local) || is.logical(local)) local <- FALSE
  local_list <- get_local_list_prediction(local)
  # interval = "prediction" needs se.fit internally to build the interval even
  # when the caller didn't request se.fit in the return value
  se_needed <- se.fit || interval == "prediction"

  model_frame <- model.frame(object)
  X <- model.matrix(object)
  y <- model.response(model_frame)
  # classical PRESS-residual identity: for ordinary (independent-error) least
  # squares, the leave-one-out residual for row i equals its ordinary
  # residual divided by (1 - leverage_i), with no need to refit n models
  cv_predict_error <- residuals(object) / (1 - hatvalues(object))
  cv_predict_val <- y - cv_predict_error

  # parallel stuff
  if (se_needed) {
    total_var <- coef(object, type = "spcov")[["ie"]]
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_se_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv_iid_se,
        vcov(object),
        Xmat = X, total_var = total_var
      )
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_se_list <- lapply(seq_len(object$n), get_loocv_iid_se, vcov(object),
        Xmat = X, total_var = total_var
      )
    }
    cv_predict_se <- vapply(cv_predict_se_list, function(x) x$se.fit, numeric(1))
  }

  bias <- mean(cv_predict_error)
  MSPE <- mean((cv_predict_error)^2)
  RMSPE <- sqrt(MSPE)
  cor2 <- cor(cv_predict_val, y)^2

  loocv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE,
    cor2 = cor2
  )

  if (interval == "prediction") {
    loocv_stats <- add_cv_coverage(loocv_stats, y, cv_predict_val, cv_predict_se, level)
  }

  if (!cv_predict && !se.fit) {
    return(loocv_stats)
  } else {
    loocv_out <- list()
    loocv_out$stats <- loocv_stats

    if (cv_predict) {
      loocv_out$cv_predict <- cv_predict_val
    }

    if (se.fit) {
      loocv_out$se.fit <- as.vector(cv_predict_se)
    }
    return(loocv_out)
  }
}

#' Append empirical prediction interval coverage to a \code{loocv()}/\code{kcv()} stats tibble
#'
#' Shared by \code{loocv.splm()}/\code{loocv.spautor()}/\code{loocv_iid()} and
#' \code{kcv.splm()}/\code{kcv.spautor()} -- only ever called when
#' \code{interval == "prediction"}, i.e. when \code{fit}/\code{se} were already
#' computed (see each caller's \code{se_needed} flag).
#'
#' @param stats_tibble The existing bias/MSPE/RMSPE/cor2 tibble
#' @param y The observed response
#' @param fit The cross-validated predictions
#' @param se The cross-validated prediction standard errors
#' @param level The prediction interval level
#'
#' @return \code{stats_tibble} with a \code{cover.XX} column appended, where
#'   \code{XX} is \code{level} with its leading \code{"0."} dropped (e.g.
#'   \code{cover.95} for \code{level = 0.95})
#'
#' @noRd
add_cv_coverage <- function(stats_tibble, y, fit, se, level) {
  # same normal-quantile interval predict.spmodel() uses by default (untouched
  # scale/df, i.e. qt(..., df = Inf) == qnorm(...))
  tstar <- qnorm(1 - (1 - level) / 2)
  lwr <- fit - tstar * se
  upr <- fit + tstar * se
  cover_name <- paste0("cover.", sub("^0\\.", "", as.character(level)))
  stats_tibble[[cover_name]] <- mean(y >= lwr & y <= upr)
  stats_tibble
}
