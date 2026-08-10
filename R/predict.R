#' Model predictions (Kriging)
#'
#' @description Predicted values and intervals based on a fitted model object.
#'
#' @param object A fitted model object.
#' @param newdata A data frame or \code{sf} object in which to
#'   look for variables with which to predict. If a data frame, \code{newdata}
#'   must contain all variables used by \code{formula(object)} and all variables
#'   representing coordinates. If an \code{sf} object, \code{newdata} must contain
#'   all variables used by \code{formula(object)} and coordinates are obtained
#'   from the geometry of \code{newdata}. If omitted, missing data from the
#'   fitted model object are used.
#' @param se.fit A logical indicating if standard errors are returned.
#'   The default is \code{FALSE}.
#' @param scale A numeric constant by which to scale the regular standard errors and intervals.
#'   Similar to but slightly different than \code{scale} for [stats::predict.lm()], because
#'   predictions form a spatial model may have different residual variances for each
#'   observation in \code{newdata}. The default is \code{NULL}, which returns
#'   the regular standard errors and intervals.
#' @param df Degrees of freedom to use for confidence or prediction intervals
#'   (ignored if \code{scale} is not specified). The default is \code{Inf}.
#' @param interval Type of interval calculation. The default is \code{"none"}.
#'   Other options are \code{"confidence"} (for confidence intervals) and
#'   \code{"prediction"} (for prediction intervals). When \code{interval}
#'   is \code{"none"} or \code{"prediction"}, predictions are returned (and when
#'   requested, their corresponding uncertainties). When \code{interval}
#'   is \code{"confidence"}, mean estimates are returned (and when
#'   requested, their corresponding uncertainties). This \code{"none"} behavior
#'   differs from that of \code{lm()}, as \code{lm()} returns confidence
#'   uncertainties (in \code{.$se.fit}).
#' @param level Tolerance/confidence level. The default is \code{0.95}.
#' @param type The prediction type, either on the response scale, link scale (only for
#'   \code{spglm()} or \code{spgautor()} model objects), terms scale,
#'   or prediction (i.e., Kriging) weight scale.
#' @param local A optional logical or list controlling the big data approximation. If omitted, \code{local}
#'   is set to \code{TRUE} or \code{FALSE} based on the observed data sample size (i.e., sample size of the fitted
#'   model object) -- if the sample size exceeds 10,000, \code{local} is
#'   set to \code{TRUE}, otherwise it is set to \code{FALSE}. This default behavior
#'   occurs because main computational
#'   burden of the big data approximation depends almost exclusively on the
#'   observed data sample size, not the number of predictions desired
#'   (which we feel is not intuitive at first glance).
#'   If \code{local} is \code{FALSE}, no big data approximation
#'   is implemented. If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{method}: The big data approximation method. If \code{method = "all"},
#'       all observations are used and \code{size} is ignored. If \code{method = "distance"},
#'       the \code{size} data observations closest (in terms of Euclidean distance)
#'       to the observation requiring prediction are used.
#'       If \code{method = "covariance"}, the \code{size} data observations
#'       with the highest covariance with the observation requiring prediction are used.
#'       If random effects and partition factors are not used in estimation and
#'       the spatial covariance function is monotone decreasing,
#'       \code{"distance"} and \code{"covariance"} are equivalent. The default
#'       is \code{"covariance"}. Only used with models fit using [splm()] or [spglm()].
#'     \item \code{size}: The number of data observations to use when \code{method}
#'       is \code{"distance"} or \code{"covariance"}. The default is 100. Only used
#'       with models fit using [splm()] or [spglm()].
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. This can significantly speed
#'       up computations even when \code{method = "all"} (i.e., no big data
#'       approximation is used), as predictions
#'       are spread out over multiple cores. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'     \item \code{byrow_threshold}: Only relevant when a random effect or partition
#'       factor is used and \code{method} is \code{"distance"} or \code{"covariance"}.
#'       In this scenario, when the observed sample size times the number of predictions
#'       is less than \code{byrow_threshold}, computations are performed once for
#'       all predictions simultaneously (which is generally faster but uses more memory).
#'       Otherwise, computations are performed separately for each prediction, one row
#'       at a time (which is generally slower but uses less memory). The default is
#'       \code{10000^2}. Only used with models fit using [splm()] or [spglm()].
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 100, method = "covariance", parallel = FALSE)}.
#'
#'   If \code{block} is \code{TRUE}, \code{local} accepts \code{method} and \code{size},
#'   and \code{method} takes values of \code{"all"}, \code{"covariance"},
#'   and \code{"distance"}, similar as when \code{block} is \code{FALSE}.
#'   The default \code{method} is \code{"covariance"} with size \code{4000}. This default \code{size} is
#'   much larger than when \code{block} is \code{FALSE}. This is because when \code{block} is \code{TRUE},
#'   covariances and explanatory variables are averaged before prediction, which greatly
#'   reduces computational burden, only requiring the Cholesky decomposition
#'   of one observed covariance matrix. Because the computational burden is reduced dramatically when \code{block} is \code{TRUE},
#'   parallel processing is not needed and hence, \code{parallel} and \code{ncores} are ignored if specified in \code{local}.
#' @param terms If \code{type} is \code{"terms"}, the type of terms to be returned,
#'   specified via either numeric position or name. The default is all terms are included.
#' @param na.action Missing (\code{NA}) values in \code{newdata} will return an error and should
#'   be removed before proceeding.
#' @param block A logical indicating whether a block prediction over the entire region
#'   in \code{newdata} should be returned. When \code{block} is \code{TRUE},
#'   \code{newdata} should be a dense grid of prediction locations that span
#'   the entire region. The default is \code{FALSE}, which
#'   returns point predictions for each location on \code{newdata}.
#' @param ... Other arguments. Only used for models fit using \code{splmRF()}
#'   or \code{spautorRF()} where \code{...} indicates other
#'   arguments to \code{ranger::predict.ranger()}.
#'
#' @details For \code{splm} and \code{spautor} objects, the (empirical) best linear unbiased predictions (i.e., Kriging
#'   predictions) at each site are returned when \code{interval} is \code{"none"}
#'   or \code{"prediction"} alongside standard errors. Prediction intervals
#'   are also returned if \code{interval} is \code{"prediction"}. When
#'   \code{interval} is \code{"confidence"}, the estimated mean is returned
#'   alongside standard errors and confidence intervals for the mean. For \code{splm_list}
#'   and \code{spautor_list} objects, predictions and associated intervals and standard errors are returned
#'   for each list element.
#'
#'   For \code{splmRF} or \code{spautorRF} objects, random forest spatial residual
#'   model predictions are computed by combining the random forest prediction with
#'   the (empirical) best linear unbiased prediction for the residual. This
#'   approach is called random forest regression Kriging. For \code{splmRF_list}
#'   or \code{spautorRF} objects,
#'   predictions are returned for each list element.
#'
#'   For \code{decorrelate} objects, the spatial decorrelation transformation
#'   predictions recorrelated to the original scale. For \code{decorrelate_list}
#'   objects, predictions are returned for each list element.
#'
#' @return For \code{splm} or \code{spautor} objects, if \code{se.fit} is \code{FALSE}, \code{predict()} returns
#'   a vector of predictions or a matrix of predictions with column names
#'   \code{fit}, \code{lwr}, and \code{upr} if \code{interval} is \code{"confidence"}
#'   or \code{"prediction"}. If \code{se.fit} is \code{TRUE}, a list with the following components is returned:
#'   \itemize{
#'     \item \code{fit}: vector or matrix as above
#'     \item \code{se.fit}: standard error of each fit
#'   }
#'
#'   For \code{splm_list} or \code{spautor_list} objects, a list that contains relevant quantities for each
#'   list element.
#'
#'   For \code{splmRF} or \code{spautorRF} objects, a vector of predictions. For \code{splmRF_list}
#'   or \code{spautorRF_list} objects, a list that contains relevant quantities for each list element.
#'
#'   For \code{decorrelate} objects, a vector of predictions. For \code{decorrelate_list}
#'   objects, a list that contains relevant quantities for each list element.
#'
#' @name predict.spmodel
#' @method predict splm
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(sulfate ~ 1, data = sulfate, spcov_type = "exponential")
#' predict(spmod, sulfate_preds)
#' predict(spmod, sulfate_preds, interval = "prediction")
#' augment(spmod, newdata = sulfate_preds, interval = "prediction")
predict.splm <- function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf, interval = c("none", "confidence", "prediction"),
                         level = 0.95, type = c("response", "terms", "weight"), block = FALSE, local, terms = NULL, na.action = na.fail, ...) {
  # match interval argument so the three display
  interval <- match.arg(interval)
  type <- match.arg(type)
  if (type == "weight") {
    se.fit <- FALSE
    interval <- "none"
  }

  # deal with local
  if (missing(local)) local <- NULL

  if (block) {
    object <- predict_block_splm(object, newdata, se.fit, scale, df, interval, level, type, local, terms, na.action, ...)
    return(object)
  }

  # build the prediction setup object and pull its elements into named local
  # variables (explicit assignment, not list2env(), so static analysis --
  # R CMD check's codetools-based check and RStudio's diagnostics -- can see
  # where obdata/xcoord/ycoord/newdata/etc. below come from)
  prediction_object <- get_prediction_object_splm(object, newdata, scale, local)
  local <- prediction_object$local
  obdata <- prediction_object$obdata
  xcoord <- prediction_object$xcoord
  ycoord <- prediction_object$ycoord
  newdata <- prediction_object$newdata
  add_newdata_rows <- prediction_object$add_newdata_rows
  spcov_params_val <- prediction_object$spcov_params_val
  randcov_params_val <- prediction_object$randcov_params_val
  newdata_model <- prediction_object$newdata_model
  offset <- prediction_object$offset

  # call terms if needed
  if (type == "terms") {
    return(predict_terms(object, newdata_model, se.fit, scale, df, interval, level, add_newdata_rows, terms, ...))
  }

  # storing newdata as a list
  npred <- NROW(newdata)
  newdata_rows_list <- split(newdata, seq_len(npred))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(npred))

  # storing newdata as a list (row_index lets get_pred_splm() slice its row
  # out of the medium-sized-data matrices built below, when applicable)
  newdata_list <- mapply(
    x = newdata_rows_list, y = newdata_model_list, i = seq_len(npred),
    FUN = function(x, y, i) list(row = x, x0 = y, row_index = i), SIMPLIFY = FALSE
  )

  if (interval %in% c("none", "prediction")) {
    # local prediction list
    local_list <- get_local_list_prediction(local)

    dotlist <- list(...)
    extra_randcov_partition_lists <- get_extra_randcov_partition_lists(object, obdata, newdata, dotlist)
    randcov_terms <- extra_randcov_partition_lists$randcov_terms
    reform_bar2 <- extra_randcov_partition_lists$reform_bar2
    partition_index_obdata <- extra_randcov_partition_lists$partition_index_obdata

    medium_precompute <- get_medium_precompute(
      object, newdata, obdata, xcoord, ycoord, npred, local_list,
      randcov_params_val, randcov_terms, spcov_params_val, reform_bar2, partition_index_obdata
    )
    dist_matrix_full <- medium_precompute$dist_matrix_full
    partition_vector_full <- medium_precompute$partition_vector_full
    cov_vector_full <- medium_precompute$cov_vector_full

    # matrix cholesky
    if (local_list$method == "all") {
      cov_matrix_val <- covmatrix(object)
      # handling closed form of none covariance
      if (inherits(spcov_params_val, c("none", "ie")) && is.null(randcov_params_val)) {
        cov_lowchol <- cov_matrix_val
        diag(cov_lowchol) <- sqrt(diag(cov_lowchol)) # already diagonal don't need transpose
      } else {
        cov_lowchol <- t(chol(cov_matrix_val))
      }
    } else {
      cov_lowchol <- NULL
    }

    # extend the prediction object with everything get_pred_splm() needs that
    # is constant across prediction rows, so it can rely on a single object
    # argument instead of two dozen individually named ones. assignment must
    # go through `[names(.)] <-` rather than `$<-`/`[[<-`,
    # since those delete a key entirely when its value is NULL (e.g. random,
    # partition_factor) instead of storing the NULL
    pred_row_context <- list(
      se.fit = se.fit, interval = interval, formula = object$terms,
      random = object$random, randcov_terms = randcov_terms,
      partition_factor = object$partition_factor, reform_bar2 = reform_bar2,
      partition_index_obdata = partition_index_obdata, cov_lowchol = cov_lowchol,
      Xmat = model.matrix(object), y = model.response(model.frame(object)),
      offset = model.offset(model.frame(object)), dim_coords = object$dim_coords,
      betahat = coefficients(object), cov_betahat = vcov(object),
      contrasts = object$contrasts, local = local_list, xlevels = object$xlevels,
      diagtol = object$diagtol, type = type,
      dist_matrix_full = dist_matrix_full, partition_vector_full = partition_vector_full,
      cov_vector_full = cov_vector_full
    )
    prediction_object[names(pred_row_context)] <- pred_row_context

    pred_splm <- run_pred_dispatch(get_pred_splm, newdata_list, prediction_object, local_list)

    if (type == "weight") {
      fit <- do.call("rbind", lapply(pred_splm, function(x) x$fit))
      if (add_newdata_rows) {
        colnames(fit) <- object$observed_index
        rownames(fit) <- object$missing_index
      }
      return(fit)
    }

    if (interval == "none") {
      fit <- vapply(pred_splm, function(x) x$fit, numeric(1))
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
      }
      if (se.fit) {
        vars <- vapply(pred_splm, function(x) x$var, numeric(1))
        se <- sqrt(vars)
        if (!is.null(scale)) {
          se <- se * scale
        }
      } else {
        se <- NULL
      }
      return(finalize_interval_none(fit, se, add_newdata_rows, object$missing_index))
    }

    if (interval == "prediction") {
      fit <- vapply(pred_splm, function(x) x$fit, numeric(1))
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
      }
      vars <- vapply(pred_splm, function(x) x$var, numeric(1))
      se <- sqrt(vars)
      if (!is.null(scale)) {
        se <- se * scale
        df <- df
      } else {
        df <- Inf
      }
      tstar <- qt(1 - (1 - level) / 2, df = df)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows, object$missing_index))
    }
  } else if (interval == "confidence") {
    # finding fitted values of the mean parameters
    fit <- as.numeric(newdata_model %*% coef(object))
    # apply offset
    if (!is.null(offset)) {
      fit <- fit + offset
    }
    newdata_model_list <- split(newdata_model, seq_len(NROW(newdata_model)))
    vars <- as.numeric(vapply(newdata_model_list, function(x) crossprod(x, vcov(object) %*% x), numeric(1)))
    se <- sqrt(vars)
    if (!is.null(scale)) {
      se <- se * scale
      df <- df
    } else {
      df <- Inf
    }
    tstar <- qt(1 - (1 - level) / 2, df = df)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows, object$missing_index))
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}

#' @rdname predict.spmodel
#' @method predict spautor
#' @order 2
#' @export
predict.spautor <- function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf, interval = c("none", "confidence", "prediction"),
                            level = 0.95, type = c("response", "terms", "weight"), local, terms = NULL, na.action = na.fail, ...) {
  # match interval argument so the three display
  interval <- match.arg(interval)
  type <- match.arg(type)
  if (type == "weight") {
    se.fit <- FALSE
    interval <- "none"
  }

  # deal with local
  if (missing(local)) local <- NULL

  # build the prediction setup object and pull its elements into named local
  # variables (explicit assignment, not list2env(), so static analysis --
  # R CMD check's codetools-based check and RStudio's diagnostics -- can see
  # where newdata/spcov_params_val/newdata_model/etc. below come from)
  prediction_object <- get_prediction_object_spautor(object, newdata, scale, local)
  local <- prediction_object$local
  newdata <- prediction_object$newdata
  spcov_params_val <- prediction_object$spcov_params_val
  randcov_params_val <- prediction_object$randcov_params_val
  newdata_model <- prediction_object$newdata_model
  offset <- prediction_object$offset

  # call terms if needed
  if (type == "terms") {
    # may want to add add_newdata_rows if we allow spautor prediction for the observed model matrix
    return(predict_terms(object, newdata_model, se.fit, scale, df, interval, level, add_newdata_rows = TRUE, terms, ...))
  }


  # storing newdata as a list
  newdata_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  if (interval %in% c("none", "prediction")) {
    # unlike splm(), spautor() models don't offer a big-data local
    # approximation here: the neighbor structure (object$W) already fixes a
    # joint covariance matrix over all sites (observed and missing), so it is
    # built once for everyone and then sliced by row/column index below,
    # rather than recomputed per prediction row
    # randcov
    randcov_Zs_val <- get_randcov_Zs(randcov_names = names(randcov_params_val), data = object$data)
    # making the partition matrix
    partition_matrix_val <- partition_matrix(object$partition_factor, object$data)
    # making the covariance matrix
    cov_matrix_val <- cov_matrix(spcov_params_val, object$W, randcov_params_val, randcov_Zs_val, partition_matrix_val, object$M)

    # making the covariance vector: the missing-by-observed block of the
    # joint covariance matrix gives each prediction row's covariance with
    # every observed row in one slice
    cov_vector_val <- cov_matrix_val[object$missing_index, object$observed_index, drop = FALSE]

    # splitting the covariance vector
    cov_vector_val_list <- split(cov_vector_val, seq_len(NROW(cov_vector_val)))

    # lower triangular cholesky
    cov_matrix_lowchol <- t(chol(cov_matrix_val[object$observed_index, object$observed_index, drop = FALSE]))

    # find X observed
    X <- model.matrix(object)
    y <- model.response(model.frame(object))
    SqrtSigInv_X <- forwardsolve(cov_matrix_lowchol, X)
    SqrtSigInv_y <- forwardsolve(cov_matrix_lowchol, y)

    # beta hat
    betahat <- coef(object)
    # residuals pearson
    residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
    # cov beta hat
    cov_betahat <- vcov(object)
    # total var
    total_var_list <- as.list(diag(cov_matrix_val[object$missing_index, object$missing_index, drop = FALSE]))


    # local prediction list (only for parallel)
    local_list <- get_local_list_prediction(local)

    # storing each new observation's row-specific quantities as a list
    cluster_list <- mapply(
      x0 = newdata_model_list, c0 = cov_vector_val_list, s0 = total_var_list,
      FUN = function(x0, c0, s0) list(x0 = x0, c0 = c0, s0 = s0), SIMPLIFY = FALSE
    )

    # extend the prediction object with everything get_pred_spautor() needs
    # that is constant across prediction rows, so it can rely on a single
    # object argument instead of many individually named ones. assignment
    # must go through `[names(.)] <-` rather than `$<-`/`[[<-`, since those
    # delete a key entirely when its value is NULL instead of storing the NULL
    pred_row_context <- list(
      cov_matrix_lowchol = cov_matrix_lowchol, betahat = betahat,
      residuals_pearson = residuals_pearson, cov_betahat = cov_betahat,
      SqrtSigInv_X = SqrtSigInv_X, se.fit = se.fit, interval = interval,
      type = type, Xmat = X
    )
    prediction_object[names(pred_row_context)] <- pred_row_context

    pred_spautor <- run_pred_dispatch(get_pred_spautor, cluster_list, prediction_object, local_list)

    if (type == "weight") {
      fit <- do.call("rbind", lapply(pred_spautor, function(x) x$fit))
      colnames(fit) <- object$observed_index
      rownames(fit) <- object$missing_index
      return(fit)
    }

    if (interval == "none") {
      fit <- vapply(pred_spautor, function(x) x$fit, numeric(1))
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
      }
      if (se.fit) {
        vars <- vapply(pred_spautor, function(x) x$var, numeric(1))
        se <- sqrt(vars)
        if (!is.null(scale)) {
          se <- se * scale
        }
      } else {
        se <- NULL
      }
      return(finalize_interval_none(fit, se, add_newdata_rows = TRUE, object$missing_index))
    }

    if (interval == "prediction") {
      fit <- vapply(pred_spautor, function(x) x$fit, numeric(1))
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
      }
      vars <- vapply(pred_spautor, function(x) x$var, numeric(1))
      se <- sqrt(vars)
      if (!is.null(scale)) {
        se <- se * scale
        df <- df
      } else {
        df <- Inf
      }
      tstar <- qt(1 - (1 - level) / 2, df = df)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows = TRUE, object$missing_index))
    }
  } else if (interval == "confidence") {
    # finding fitted values of the mean parameters
    fit <- as.numeric(newdata_model %*% coef(object))
    # apply offset
    if (!is.null(offset)) {
      fit <- fit + offset
    }
    vars <- as.numeric(vapply(newdata_model_list, function(x) crossprod(x, vcov(object) %*% x), numeric(1)))
    se <- sqrt(vars)
    if (!is.null(scale)) {
      se <- se * scale
      df <- df
    } else {
      df <- Inf
    }
    tstar <- qt(1 - (1 - level) / 2, df = df)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows = TRUE, object$missing_index))
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}

#' Predict a single new observation for an \code{splm()} model
#'
#' @param newdata_list A list with elements \code{row} (a single-row data frame
#'   for the new observation), \code{x0} (its design matrix row), and
#'   \code{row_index} (its row number within \code{newdata}, used to slice
#'   \code{dist_matrix_full}/\code{partition_vector_full}/\code{cov_vector_full}
#'   when those are supplied)
#' @param prediction_object A list of values constant across prediction rows
#'   (built by \code{predict.splm()}), assigned into named local variables at
#'   the top of this function: \code{se.fit}, \code{interval}, \code{formula},
#'   \code{obdata} (the observed data, or, for the big-data local methods, the
#'   full observed data before neighbor subsetting), \code{xcoord}, \code{ycoord},
#'   \code{spcov_params_val}, \code{random} (a random effect formula, or
#'   \code{NULL}), \code{randcov_params_val} (a \code{randcov_params} object,
#'   or \code{NULL}), \code{randcov_terms} (see \code{get_extra_randcov_list()}),
#'   \code{partition_factor} (a partition factor formula, or \code{NULL}),
#'   \code{reform_bar2} (the partition factor's grouping formula, or
#'   \code{NULL}), \code{partition_index_obdata} (see
#'   \code{get_extra_partition_list()}), \code{cov_lowchol} (the lower Cholesky
#'   factor of the observed-data covariance matrix, ignored and recomputed
#'   from the local neighborhood when \code{local$method} is \code{"distance"}
#'   or \code{"covariance"}), \code{Xmat} (the observed-data design matrix),
#'   \code{y} (the observed response vector), \code{offset} (the observed-data
#'   offset, or \code{NULL}), \code{betahat}, \code{cov_betahat}, \code{dim_coords},
#'   \code{contrasts}, \code{local} (a fully-specified big-data \code{local}
#'   list), \code{xlevels}, \code{diagtol}, \code{type} (\code{"response"}
#'   or \code{"weight"}), and, for medium-sized data (see \code{predict.splm()}),
#'   \code{dist_matrix_full}, \code{partition_vector_full}, and
#'   \code{cov_vector_full} -- the observed-by-prediction distance, partition,
#'   and covariance matrices, precomputed once (vectorized across every
#'   prediction row) so this row's values can be sliced out of them instead of
#'   reconstructed from scratch; \code{NULL} for big data, where reconstructing
#'   them one row at a time instead of all at once avoids ever holding an
#'   observed-by-prediction matrix in memory
#'
#' @return A list with element \code{fit} (and, if \code{se.fit} or
#'   \code{interval == "prediction"}, element \code{var}) for the new observation
#'
#' @noRd
get_pred_splm <- function(newdata_list, prediction_object) {
  # explicit assignment, not list2env(), so static analysis (R CMD check's
  # codetools-based check and RStudio's diagnostics) can see where each name
  # below comes from
  se.fit <- prediction_object$se.fit
  interval <- prediction_object$interval
  formula <- prediction_object$formula
  obdata <- prediction_object$obdata
  xcoord <- prediction_object$xcoord
  ycoord <- prediction_object$ycoord
  spcov_params_val <- prediction_object$spcov_params_val
  random <- prediction_object$random
  randcov_params_val <- prediction_object$randcov_params_val
  randcov_terms <- prediction_object$randcov_terms
  partition_factor <- prediction_object$partition_factor
  reform_bar2 <- prediction_object$reform_bar2
  partition_index_obdata <- prediction_object$partition_index_obdata
  cov_lowchol <- prediction_object$cov_lowchol
  Xmat <- prediction_object$Xmat
  y <- prediction_object$y
  offset <- prediction_object$offset
  betahat <- prediction_object$betahat
  cov_betahat <- prediction_object$cov_betahat
  dim_coords <- prediction_object$dim_coords
  contrasts <- prediction_object$contrasts
  local <- prediction_object$local
  xlevels <- prediction_object$xlevels
  diagtol <- prediction_object$diagtol
  type <- prediction_object$type
  dist_matrix_full <- prediction_object$dist_matrix_full
  partition_vector_full <- prediction_object$partition_vector_full
  cov_vector_full <- prediction_object$cov_vector_full

  # medium-mode reuse, partition-index subsetting, and dense-mode recompute
  # of the local distance/covariance vector -- see get_pred_local_setup() in
  # R/predict_helpers.R
  local_setup <- get_pred_local_setup(
    newdata_list, obdata, xcoord, ycoord, dim_coords,
    spcov_params_val, randcov_params_val, randcov_terms,
    partition_factor, reform_bar2, partition_index_obdata,
    random, local, dist_matrix_full, partition_vector_full, cov_vector_full
  )
  obdata <- local_setup$obdata
  randcov_terms <- local_setup$randcov_terms
  dist_vector <- local_setup$dist_vector
  cov_vector_val <- local_setup$cov_vector_val

  # subsetting data if method distance
  if (local$method == "distance") {
    n <- length(cov_vector_val)
    # want the smallest distance here and order goes from smallest first to largest last (keep last values with are smallest distance)
    nn_index <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, local$size))]
    obdata <- obdata[nn_index, , drop = FALSE]
    cov_vector_val <- cov_vector_val[nn_index]
  }

  if (local$method == "covariance") {
    n <- length(cov_vector_val)
    # want the largest covariance here and order goes from smallest first to largest last (keep last values which are largest covariance)
    # use abs() here for largest absolute covariance
    # generally the same as cov unless spcov type 
    # is not a monotonic function of distance
    cov_index <- order(abs(as.numeric(cov_vector_val)))[seq(from = n, to = max(1, n - local$size + 1))]
    obdata <- obdata[cov_index, , drop = FALSE]
    cov_vector_val <- cov_vector_val[cov_index]
  }

  if (local$method %in% c("distance", "covariance")) {
    # this is the observed-by-observed covariance among just the retained
    # local neighborhood, which is never part of the medium-mode precomputed
    # matrices (those are observed-by-prediction only), so it must always be
    # built fresh here regardless of medium_mode
    xlev_list <- if (is.null(random)) NULL else lapply(randcov_terms, function(x) x$xlev)
    dist_matrix <- spdist(obdata, xcoord, ycoord, sparse = FALSE)
    cov_matrix_val <- get_obs_cov_matrix(
      dist_matrix, obdata, spcov_params_val, randcov_params_val,
      random, partition_factor,
      diagtol = diagtol, xlev_list = xlev_list
    )
    # the local neighborhood is small (local$size) and effectively dense (a
    # spatial covariance has no exact zeros), so factoring it as a plain base
    # matrix instead of a sparse Matrix-class object avoids paying repeated
    # S4 dispatch/validity-check overhead for a matrix with no sparsity to
    # exploit -- this branch runs once per prediction row, so that overhead
    # is what dominated profiling at scale
    cov_lowchol <- t(base::chol(as.matrix(cov_matrix_val)))
    model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass, xlev = xlevels)
    Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
    y <- model.response(model_frame)
    offset <- model.offset(model_frame)
  }


  # handle offset
  if (!is.null(offset)) {
    y <- y - offset
  }

  # "whiten" X, y, and c0 by left-multiplying by Sigma^{-1/2} -- implemented
  # as forward substitution against the lower Cholesky factor rather than
  # explicitly forming Sigma^{-1}, which is both cheaper and more numerically
  # stable
  c0 <- as.numeric(cov_vector_val)
  SqrtSigInv_X <- base::forwardsolve(cov_lowchol, Xmat)
  SqrtSigInv_y <- base::forwardsolve(cov_lowchol, y)
  SqrtSigInv_c0 <- base::forwardsolve(cov_lowchol, c0)
  x0 <- newdata_list$x0

  if (type == "weight") {
    Xt_SigInv <- base::t(base::backsolve(base::t(cov_lowchol), SqrtSigInv_X))
    betahat_wt <- cov_betahat %*% Xt_SigInv # for big data, for this to exactly equal wts %*% y = "fit" for
    # type != weight, Xt_SigInv should use SigInv from the original fit, which would require recomputing cholprods.
    # For now, these are approximate.
    residuals_weight <- -1 * Xmat %*% betahat_wt # this is recomputed over and over when using all data consider making more efficient
    diag(residuals_weight) <- diag(residuals_weight) + 1
    fit <- x0 %*% betahat_wt + base::crossprod(SqrtSigInv_c0, base::forwardsolve(cov_lowchol, residuals_weight))
    if (local$method %in% c("distance", "covariance")) {
      wtfit <- fit
      fit <- Matrix::Matrix(0, nrow = 1, ncol = n, sparse = TRUE)
      if (local$method == "distance") {
        fit[nn_index] <- wtfit
      }
      if (local$method == "covariance") {
        fit[cov_index] <- wtfit
      }
    }
  } else {
    # universal kriging BLUP: the trend x0 %*% betahat plus a covariance-
    # weighted combination of the observed (whitened) residuals -- the
    # closer/more correlated an observation is with the new location (larger
    # entries of c0), the more its residual pulls the prediction away from
    # the trend line
    residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
    fit <- as.numeric(x0 %*% betahat + base::crossprod(SqrtSigInv_c0, residuals_pearson))
  }
  if (se.fit || interval == "prediction") {
    # kriging prediction variance = marginal variance of the new observation,
    # minus the variance explained by conditioning on the observed data
    # (crossprod(SqrtSigInv_c0, SqrtSigInv_c0)), plus an inflation term
    # (H %*% cov_betahat %*% t(H)) accounting for betahat itself being
    # estimated rather than known
    H <- x0 - base::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
    # a random slope's contribution to Var(Y0) is sigma^2 * x0^2, not sigma^2
    # (as it would be for a random intercept), so it must be computed for this
    # specific newdata row rather than summed directly from randcov_params_val
    total_var <- spcov_params_val[["de"]] + spcov_params_val[["ie"]] +
      randcov_newvar(randcov_params_val, newdata_list$row, randcov_terms)
    var <- as.numeric(total_var - base::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% base::tcrossprod(cov_betahat, H))
    pred_list <- list(fit = fit, var = var)
  } else {
    pred_list <- list(fit = fit)
  }
  pred_list
}

#' Predict a single new observation for an \code{spautor()} model
#'
#' @param cluster_list A list with elements \code{x0} (the design matrix row
#'   for the new observation), \code{c0} (the covariance vector between the
#'   new observation and the observed data), and \code{s0} (the marginal
#'   variance of the new observation)
#' @param prediction_object A list of values constant across prediction rows
#'   (built by \code{predict.spautor()}), assigned into named local variables
#'   at the top of this function: \code{cov_matrix_lowchol} (the lower Cholesky factor of
#'   the observed-data covariance matrix), \code{betahat}, \code{residuals_pearson}
#'   (Pearson residuals for the observed data), \code{cov_betahat}, \code{SqrtSigInv_X}
#'   (\code{cov_matrix_lowchol^{-1} \%*\% Xmat}), \code{se.fit}, \code{interval}
#'   (\code{"none"}, \code{"confidence"}, or \code{"prediction"}), \code{type}
#'   (\code{"response"} or \code{"weight"}), and \code{Xmat} (the observed-data
#'   design matrix)
#'
#' @return A list with element \code{fit} (and, if \code{se.fit} or
#'   \code{interval == "prediction"}, element \code{var}) for the new observation
#'
#' @noRd
get_pred_spautor <- function(cluster_list, prediction_object) {
  # explicit assignment, not list2env(), so static analysis (R CMD check's
  # codetools-based check and RStudio's diagnostics) can see where each name
  # below comes from
  cov_matrix_lowchol <- prediction_object$cov_matrix_lowchol
  betahat <- prediction_object$betahat
  residuals_pearson <- prediction_object$residuals_pearson
  cov_betahat <- prediction_object$cov_betahat
  SqrtSigInv_X <- prediction_object$SqrtSigInv_X
  se.fit <- prediction_object$se.fit
  interval <- prediction_object$interval
  type <- prediction_object$type
  Xmat <- prediction_object$Xmat

  # shared with get_pred_spgautor() -- see get_areal_pred() in
  # R/predict_helpers.R
  get_areal_pred(cluster_list, cov_matrix_lowchol, betahat, residuals_pearson, cov_betahat, SqrtSigInv_X, se.fit, interval, type, Xmat)
}

#' @name predict.spmodel
#' @method predict splm_list
#' @order 3
#' @export
predict.splm_list <- function(object, newdata, se.fit = FALSE, scale = NULL,
                              df = Inf, interval = c("none", "confidence", "prediction"),
                              level = 0.95, type = c("response", "terms"), local,
                              terms = NULL, na.action = na.fail, ...) {
  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with local
  if (missing(local)) local <- NULL

  if (missing(newdata)) {
    preds <- lapply(object, function(x) {
      predict(
        x,
        se.fit = se.fit,
        scale = scale,
        df = df,
        interval = interval,
        level = level,
        type = type,
        local = local,
        terms = terms, # don't need na.action as it is fixed in predict
        ...
      )
    })
  } else {
    preds <- lapply(object, function(x) {
      predict(
        x,
        newdata = newdata,
        se.fit = se.fit,
        scale = scale,
        df = df,
        interval = interval,
        level = level,
        type = type,
        local = local,
        terms = terms, # don't need na.action as it is fixed in predict
        ...
      )
    })
  }
  names(preds) <- names(object)
  preds
}

#' @name predict.spmodel
#' @method predict spautor_list
#' @order 4
#' @export
predict.spautor_list <- predict.splm_list


#' @rdname predict.spmodel
#' @method predict splmRF
#' @order 5
#' @export
#'
#' @references
#' Fox, E.W., Ver Hoef, J. M., & Olsen, A. R. (2020). Comparing spatial
#'   regression to random forests for large environmental data sets.
#'   \emph{PloS one}, 15(3), e0229509.
#'
#' @examples
#' \donttest{
#' sulfate$var <- rnorm(NROW(sulfate)) # add noise variable
#' sulfate_preds$var <- rnorm(NROW(sulfate_preds)) # add noise variable
#' sprfmod <- splmRF(sulfate ~ var, data = sulfate, spcov_type = "exponential")
#' predict(sprfmod, sulfate_preds)
#' }
predict.splmRF <- function(object, newdata, local, ...) {
  # check to see if ranger installed
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using predict with an splmRF or spautorRF object", call. = FALSE)
  } else {
    # find newdata if required
    if ((missing(newdata) && !is.null(object$newdata))) {
      newdata <- object$newdata
    }

    # get ... objects
    dotlist <- as.list(substitute(alist(...)))[-1]
    dotlist_names <- names(dotlist)

    # hardcode ranger names because of predict.ranger export issue
    ranger_names <- c("predict.all", "num.trees", "se.method", "quantiles", "what", "seed", "num.threads", "verbose")
    ranger_args <- dotlist[dotlist_names %in% ranger_names]

    # do random forest prediction
    ranger_pred <- do.call(predict, c(list(object = object$ranger, data = as.data.frame(newdata), type = "response"), ranger_args))

    # set local if missing
    if (missing(local)) local <- NULL
    # do splm prediction
    splm_pred <- do.call(predict, list(object = object$splm, newdata = newdata, local = local))
  }
  # obtain final predictions
  ranger_pred$predictions + splm_pred
}

#' @rdname predict.spmodel
#' @method predict spautorRF
#' @order 6
#' @export
predict.spautorRF <- function(object, newdata, local, ...) {
  # check to see if ranger installed
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using predict with an splmRF or spautorRF object", call. = FALSE)
  } else {
    # find newdata
    newdata <- object$newdata


    # get ... objects
    dotlist <- as.list(substitute(alist(...)))[-1]
    dotlist_names <- names(dotlist)

    # hardcode ranger names because of predict.ranger export issue
    ranger_names <- c("predict.all", "num.trees", "se.method", "quantiles", "what", "seed", "num.threads", "verbose")
    ranger_args <- dotlist[dotlist_names %in% ranger_names]

    # do random forest prediction
    ranger_pred <- do.call(predict, c(list(object = object$ranger, data = as.data.frame(newdata), type = "response"), ranger_args))

    # set local if missing
    if (missing(local)) local <- NULL
    # do spautor prediction -- newdata is intentionally omitted here (unlike
    # the splmRF/ranger call above): object$spautor already knows its own
    # prediction locations via object$spautor$newdata, which predict.spautor()
    # requires newdata (if supplied at all) to match exactly, and this
    # wrapper's own newdata is not identical to it (e.g. sf geometry handling
    # differs) even though both represent the same underlying rows
    spautor_pred <- do.call(predict, list(object = object$spautor, local = local))
  }
  # obtain final predictions
  ranger_pred$predictions + spautor_pred
}

#' @name predict.spmodel
#' @method predict splmRF_list
#' @order 7
#' @export
predict.splmRF_list <- function(object, newdata, local, ...) {
  # check to see if ranger installed
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using predict with an splmRF_list object", call. = FALSE)
  } else {
    # find newdata if required
    if ((missing(newdata) && !is.null(object$newdata))) {
      newdata <- object$newdata
    }

    # get ... objects
    dotlist <- as.list(substitute(alist(...)))[-1]
    dotlist_names <- names(dotlist)

    # hardcode ranger names because of predict.ranger export issue
    ranger_names <- c("predict.all", "num.trees", "se.method", "quantiles", "what", "seed", "num.threads", "verbose")
    ranger_args <- dotlist[dotlist_names %in% ranger_names]

    # do random forest prediction
    ranger_pred <- do.call(predict, c(list(object = object$ranger, data = as.data.frame(newdata), type = "response"), ranger_args))

    # set local if missing
    if (missing(local)) local <- NULL
    # do splm prediction
    splm_pred <- lapply(object$splm_list, function(x) {
      do.call(predict, list(object = x, newdata = newdata, local = local))
    })
  }
  # obtain final predictions
  sprf_pred <- lapply(splm_pred, function(x) ranger_pred$predictions + x)
  names(sprf_pred) <- names(object$splm_list)
  sprf_pred
}

#' @name predict.spmodel
#' @method predict spautorRF_list
#' @order 8
#' @export
predict.spautorRF_list <- function(object, newdata, local, ...) {
  # check to see if ranger installed
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using predict with an splmRF_list object", call. = FALSE)
  } else {
    # find newdata
    newdata <- object$newdata

    # get ... objects
    dotlist <- as.list(substitute(alist(...)))[-1]
    dotlist_names <- names(dotlist)

    # hardcode ranger names because of predict.ranger export issue
    ranger_names <- c("predict.all", "num.trees", "se.method", "quantiles", "what", "seed", "num.threads", "verbose")
    ranger_args <- dotlist[dotlist_names %in% ranger_names]

    # do random forest prediction
    ranger_pred <- do.call(predict, c(list(object = object$ranger, data = as.data.frame(newdata), type = "response"), ranger_args))

    # set local if missing
    if (missing(local)) local <- NULL
    # do spautor prediction -- newdata intentionally omitted; see the matching
    # note in predict.spautorRF()
    spautor_pred <- lapply(object$spautor_list, function(x) {
      do.call(predict, list(object = x, local = local))
    })
  }
  # obtain final predictions
  sprf_pred <- lapply(spautor_pred, function(x) ranger_pred$predictions + x)
  names(sprf_pred) <- names(object$spautor_list)
  sprf_pred
}


#' Precompute per-term random effect context for prediction/loocv
#'
#' @param object A fitted model object from [splm()] or [spglm()]
#' @param obdata The observed data
#' @param newdata The data requiring prediction
#'
#' @return A list with element \code{randcov_terms}: \code{NULL} if there are
#'   no random effects, otherwise a named list (one element per random effect
#'   term) with \code{reform_bar2}, \code{reform_bar1}, \code{group_label},
#'   \code{xlev}, \code{level_index_map}, and \code{slope_val}, computed once
#'   per \code{predict()}/\code{loocv()} call so downstream code (e.g.
#'   \code{get_randcov_vectors()}) can look up matches in O(1) instead of
#'   rescanning \code{obdata} for every prediction row
#'
#' @noRd
get_extra_randcov_list <- function(object, obdata, newdata) {
  # random stuff: one context object per random effect term (reform_bar2,
  # reform_bar1, group_label, xlev, level_index_map, slope_val) instead of
  # four separate same-keyed lists threaded through every downstream call
  if (is.null(object$random)) {
    return(list(randcov_terms = NULL))
  }
  randcov_names <- get_randcov_names(object$random)
  randcov_terms <- lapply(randcov_names, function(randcov_name) {
    bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)

    reform_bar2_mf <- model.frame(reform_bar2, obdata)
    reform_bar2_terms <- terms(reform_bar2_mf)
    xlev <- .getXlevels(reform_bar2_terms, reform_bar2_mf)
    group_label <- model_matrix_group_labels(reform_bar2, obdata)

    # adding dummy levels if newdata observations of random effects are not in original data
    # terms object is unchanged if levels change
    xlev_full <- .getXlevels(reform_bar2_terms, rbind(reform_bar2_mf, model.frame(reform_bar2, newdata)))
    if (!identical(xlev, xlev_full)) {
      xlev <- xlev_full
    }

    # precomputed once per predict() call (not once per prediction row): maps
    # each observed group label to the obdata row indices that carry it, so
    # get_randcov_vectors() can look up matches in O(1) instead of scanning
    # all of obdata for every prediction row
    level_index_map <- split(seq_along(group_label), group_label)

    if (bar_split[[1]] != "1") {
      reform_bar1 <- reformulate(bar_split[[1]], intercept = FALSE)
      slope_val <- as.vector(model.matrix(reform_bar1, obdata))
    } else {
      reform_bar1 <- NULL
      slope_val <- NULL
    }

    list(
      reform_bar2 = reform_bar2, reform_bar1 = reform_bar1, group_label = group_label,
      xlev = xlev, level_index_map = level_index_map, slope_val = slope_val
    )
  })
  names(randcov_terms) <- randcov_names
  list(randcov_terms = randcov_terms)
}

#' Precompute partition factor context for prediction/loocv
#'
#' @param object A fitted model object from [splm()] or [spglm()]
#' @param obdata The observed data
#' @param newdata The data requiring prediction
#'
#' @return A list with elements \code{reform_bar2} (the partition factor's
#'   grouping formula, or \code{NULL} if there is no partition factor) and
#'   \code{partition_index_obdata} (a list with \code{group_label}, \code{xlev},
#'   and \code{level_index_map} for \code{obdata}, mirroring the random effect
#'   context built by \code{get_extra_randcov_list()})
#'
#' @noRd
get_extra_partition_list <- function(object, obdata, newdata) {
  # partition factor context: group_label, xlev, and a precomputed
  # level_index_map (mirrors the random effect context in
  # get_extra_randcov_list(), built once per predict() call rather than
  # rescanning all of obdata for every prediction row)
  if (is.null(object$partition_factor)) {
    return(list(reform_bar2 = NULL, partition_index_obdata = NULL))
  }
  partition_factor_val <- get_partition_name(labels(terms(object$partition_factor)))
  bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
  reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)

  p_reform_bar2_mf <- model.frame(reform_bar2, obdata)
  p_reform_bar2_terms <- terms(p_reform_bar2_mf)
  xlev <- .getXlevels(p_reform_bar2_terms, p_reform_bar2_mf)
  group_label <- model_matrix_group_labels(reform_bar2, obdata)

  # adding dummy levels if newdata observations of the partition factor are not in original data
  # terms object is unchanged if levels change
  xlev_full <- .getXlevels(p_reform_bar2_terms, rbind(p_reform_bar2_mf, model.frame(reform_bar2, newdata)))
  if (!identical(xlev, xlev_full)) {
    xlev <- xlev_full
  }

  level_index_map <- split(seq_along(group_label), group_label)

  partition_index_obdata <- list(group_label = group_label, xlev = xlev, level_index_map = level_index_map)
  list(reform_bar2 = reform_bar2, partition_index_obdata = partition_index_obdata)
}

#' Replace unseen factor/character levels in \code{newdata} with a placeholder
#'
#' @param varnames Names of variables to check
#' @param obdata The observed data
#' @param newdata The data requiring prediction
#'
#' @return \code{newdata} with any factor or character values in \code{varnames}
#'   not present in \code{obdata} replaced by a placeholder level, so that
#'   downstream code building a design matrix does not error on an unseen
#'   level (numeric/integer variables, e.g. a continuous random slope, are left
#'   unchanged since there is no notion of an "unseen level" for them)
#'
#' @details \code{NA %in% x} is always \code{FALSE}, so without an explicit
#'   check an \code{NA} value would be silently treated as just another
#'   unseen level and replaced by the placeholder. This gives that row a zero
#'   random effect/partition factor contribution with no error or warning,
#'   rather than surfacing a "Cannot have NA values in predictors."
#'   error already used for \code{NA} in a fixed effect or random slope value.
#'
#' @noRd
replace_newdata <- function(varnames, obdata, newdata) {
  newdata_vec <- lapply(varnames, function(x) {
    obdata_vec <- obdata[[x]]
    newdata_vec <- newdata[[x]]
    if (anyNA(newdata_vec)) {
      stop("Cannot have NA values in predictors.", call. = FALSE)
    }
    index_vec <- !newdata_vec %in% obdata_vec
    if (any(index_vec)) {
      if (is.factor(newdata_vec)) {
        newdata_vec <- as.character(newdata_vec)
        newdata_vec[index_vec] <- "...this_is_a_new_level..."
        newdata_vec <- as.factor(newdata_vec)
      } else if (is.character(newdata_vec)) {
        newdata_vec[index_vec] <- "...this_is_a_new_level..."
      } else { # don't replace if numeric/integer (continuous random slope)
        newdata_vec <- NULL
      }
    }
    newdata_vec
  })
  names(newdata_vec) <- varnames
  for (x in varnames) {
    if (!is.null(newdata_vec[[x]])) {
      newdata[[x]] <- newdata_vec[[x]]
    }
  }
  newdata
}
