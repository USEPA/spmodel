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
#'   requested, their corresponding uncertainties).
#' @param level Tolerance/confidence level. The default is \code{0.95}.
#' @param type The prediction type, either on the response scale, link scale (only for
#'   \code{spglm()} or \code{spgautor()} model objects), or terms scale.
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
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 100, method = "covariance", parallel = FALSE)}.
#' @param terms If \code{type} is \code{"terms"}, the type of terms to be returned,
#'   specified via either numeric position or name. The default is all terms are included.
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
#'   the (empirical) best linear unbiased prediction for the residual. Fox et al. (2020)
#'   call this approach random forest regression Kriging. For \code{splmRF_list}
#'   or \code{spautorRF} objects,
#'   predictions are returned for each list element.
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
#' @name predict.spmodel
#' @method predict splm
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(sulfate ~ 1,
#'   data = sulfate,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' predict(spmod, sulfate_preds)
#' predict(spmod, sulfate_preds, interval = "prediction")
#' augment(spmod, newdata = sulfate_preds, interval = "prediction")
predict.splm <- function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf, interval = c("none", "confidence", "prediction"),
                         level = 0.95, type = c("response", "terms"), local, terms = NULL, ...) {

  # match interval argument so the three display
  interval <- match.arg(interval)
  type <- match.arg(type)

  # deal with local
  if (missing(local)) {
    local <- NULL
  }

  # check scale is numeric (if specified)
  if (!is.null(scale) && !is.numeric(scale)) {
    stop("scale must be numeric.", call. = FALSE)
  }

  # handle na action -- this is an inefficient workaround that should be fixed later
  # placeholder as a reminder to consider adding na.action argument at a later date
  # na_action <- as.character(substitute(na.action))

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }

  # rename relevant quantities
  obdata <- object$obdata
  xcoord <- object$xcoord
  ycoord <- object$ycoord

  # write newdata if predicting missing data
  if (missing(newdata)) {
    add_newdata_rows <- TRUE
    newdata <- object$newdata
  } else {
    add_newdata_rows <- FALSE
  }

  # deal with local
  if (is.null(local)) {
    if (object$n > 10000) {
    # if (object$n > 5000 || NROW(newdata) > 5000) {
      local <- TRUE
      message("Because the sample size of the fitted model object exceeds 10,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun predict() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save randcov param vector
  randcov_params_val <- coef(object, type = "randcov")


  attr_sp <- attr(class(newdata), "package")
  if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
    stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
  }

  if (inherits(newdata, "sf")) {
    newdata <- suppressWarnings(sf::st_centroid(newdata))

    newdata <- sf_to_df(newdata)
    names(newdata)[[which(names(newdata) == ".xcoord")]] <- as.character(xcoord) # only relevant if newdata is sf data is not
    names(newdata)[[which(names(newdata) == ".ycoord")]] <- as.character(ycoord) # only relevant if newdata is sf data is not
  }

  # add back in zero column to cover anisotropy (should make anisotropy only available 1-d)
  if (object$dim_coords == 1) {
    obdata[[ycoord]] <- 0
    newdata[[ycoord]] <- 0
  }

  if (object$anisotropy) { # could just do rotate != 0 || scale != 1
    obdata_aniscoords <- transform_anis(obdata, xcoord, ycoord,
      rotate = spcov_params_val[["rotate"]],
      scale = spcov_params_val[["scale"]]
    )
    obdata[[xcoord]] <- obdata_aniscoords$xcoord_val
    obdata[[ycoord]] <- obdata_aniscoords$ycoord_val
    newdata_aniscoords <- transform_anis(newdata, xcoord, ycoord,
      rotate = spcov_params_val[["rotate"]],
      scale = spcov_params_val[["scale"]]
    )
    newdata[[xcoord]] <- newdata_aniscoords$xcoord_val
    newdata[[ycoord]] <- newdata_aniscoords$ycoord_val
  }

  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  if (any(grepl("nmatrix.", attributes(formula_newdata)$dataClasses, fixed = TRUE)) && NROW(newdata) == 1) {
    newdata <- newdata[c(1, 1), , drop = FALSE]
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
    newdata_model <- newdata_model[1, , drop = FALSE]
    # find offset
    offset <- model.offset(newdata_model_frame)
    if (!is.null(offset)) {
      offset <- offset[1]
    }
    newdata <- newdata[1, , drop = FALSE]
  } else {
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    # assumes that predicted observations are not outside the factor levels
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
    # find offset
    offset <- model.offset(newdata_model_frame)
  }
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # call terms if needed
  if (type == "terms") {
    return(predict_terms(object, newdata_model, se.fit, scale, df, interval, level, add_newdata_rows, terms, ...))
  }

  # storing newdata as a list
  newdata_rows_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_list <- mapply(x = newdata_rows_list, y = newdata_model_list, FUN = function(x, y) list(row = x, x0 = y), SIMPLIFY = FALSE)

  if (interval %in% c("none", "prediction")) {

    # local prediction list
    local_list <- get_local_list_prediction(local)

    # hardcoding none covariance (commented out 0.5.1 for efficiency's sake in use with local)
    # if (inherits(spcov_params_val, "none") && is.null(randcov_params_val)) {
    #   local_list$method <- "all"
    # }

    dotlist <- list(...)
    dotlist_names <- names(dotlist)

    if ("extra_randcov_list" %in% dotlist_names && !is.null(dotlist[["extra_randcov_list"]])) {
      extra_randcov_list <- dotlist$extra_randcov_list
    } else {
      extra_randcov_list <- get_extra_randcov_list(object, obdata, newdata)
    }
    reform_bar2_list <- extra_randcov_list$reform_bar2_list
    Z_index_obdata_list <- extra_randcov_list$Z_index_obdata_list
    reform_bar1_list <- extra_randcov_list$reform_bar1_list
    Z_val_obdata_list <- extra_randcov_list$Z_val_obdata_list

    if ("extra_partition_list" %in% dotlist_names && !is.null(dotlist[["extra_partition_list"]])) {
      extra_partition_list <- dotlist$extra_partition_list
    } else {
      extra_partition_list <- get_extra_partition_list(object, obdata, newdata)
    }
    reform_bar2 <- extra_partition_list$reform_bar2
    partition_index_obdata <- extra_partition_list$partition_index_obdata

    # # random stuff
    # if (!is.null(object$random)) {
    #   randcov_names <- get_randcov_names(object$random)
    #   # this causes a memory leak and was not even needed
    #   # randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
    #   # comment out here for simple
    #   reform_bar_list <- lapply(randcov_names, function(randcov_name) {
    #     bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
    #     reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
    #     if (bar_split[[1]] != "1") {
    #       reform_bar1 <- reformulate(bar_split[[1]], intercept = FALSE)
    #     } else {
    #       reform_bar1 <- NULL
    #     }
    #     list(reform_bar2 = reform_bar2, reform_bar1 = reform_bar1)
    #   })
    #   reform_bar2_list <- lapply(reform_bar_list, function(x) x$reform_bar2)
    #   names(reform_bar2_list) <- randcov_names
    #   reform_bar1_list <- lapply(reform_bar_list, function(x) x$reform_bar1)
    #   names(reform_bar1_list) <- randcov_names
    #   Z_index_obdata_list <- lapply(reform_bar2_list, function(reform_bar2) {
    #
    #     reform_bar2_mf <- model.frame(reform_bar2, obdata)
    #     reform_bar2_terms <- terms(reform_bar2_mf)
    #     reform_bar2_xlev <- .getXlevels(reform_bar2_terms, reform_bar2_mf)
    #     reform_bar2_mx <- model.matrix(reform_bar2, obdata)
    #     reform_bar2_names <- colnames(reform_bar2_mx)
    #     reform_bar2_split <- split(reform_bar2_mx, seq_len(NROW(reform_bar2_mx)))
    #     reform_bar2_vals <- reform_bar2_names[vapply(reform_bar2_split, function(y) which(as.logical(y)), numeric(1))]
    #
    #     # adding dummy levels if newdata observations of random effects are not in original data
    #     # terms object is unchanged if levels change
    #     # reform_bar2_mf_new <- model.frame(reform_bar2, newdata)
    #     # reform_bar2_mf_full <- model.frame(reform_bar2, merge(obdata, newdata, all = TRUE))
    #     # reform_bar2_terms_full <- terms(rbind(reform_bar2_mf, reform_bar2_mf_new))
    #     reform_bar2_xlev_full <- .getXlevels(reform_bar2_terms, rbind(reform_bar2_mf, model.frame(reform_bar2, newdata)))
    #     if (!identical(reform_bar2_xlev, reform_bar2_xlev_full)) {
    #       reform_bar2_xlev <- reform_bar2_xlev_full
    #     }
    #
    #     list(reform_bar2_vals = reform_bar2_vals, reform_bar2_xlev = reform_bar2_xlev)
    #   })
    #   # Z_index_obdata_list <- lapply(reform_bar2_list, function(reform_bar2) as.vector(model.matrix(reform_bar2, obdata)))
    #   names(Z_index_obdata_list) <- randcov_names
    #   Z_val_obdata_list <- lapply(reform_bar1_list, function(reform_bar1) {
    #     if (is.null(reform_bar1)) {
    #       return(NULL)
    #     } else {
    #       return(as.vector(model.matrix(reform_bar1, obdata)))
    #     }
    #   })
    #   names(Z_val_obdata_list) <- randcov_names
    # } else {
    #   reform_bar2_list <- NULL
    #   Z_index_obdata_list <- NULL
    #   reform_bar1_list <- NULL
    #   Z_val_obdata_list <- NULL
    # }
    #
    # # partition factor stuff
    # if (!is.null(object$partition_factor)) {
    #   partition_factor_val <- get_partition_name(labels(terms(object$partition_factor)))
    #   bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
    #   reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
    #   p_reform_bar2_mf <- model.frame(reform_bar2, obdata)
    #   p_reform_bar2_terms <- terms(p_reform_bar2_mf)
    #   p_reform_bar2_xlev <- .getXlevels(p_reform_bar2_terms, p_reform_bar2_mf)
    #   p_reform_bar2_mx <- model.matrix(reform_bar2, obdata)
    #   p_reform_bar2_names <- colnames(p_reform_bar2_mx)
    #   p_reform_bar2_split <- split(p_reform_bar2_mx, seq_len(NROW(p_reform_bar2_mx)))
    #   p_reform_bar2_vals <- p_reform_bar2_names[vapply(p_reform_bar2_split, function(y) which(as.logical(y)), numeric(1))]
    #
    #
    #   # adding dummy levels if newdata observations of random effects are not in original data
    #   # terms object is unchanged if levels change
    #   # p_reform_bar2_mf_new <- model.frame(reform_bar2, newdata)
    #   # reform_bar2_mf_full <- model.frame(reform_bar2, merge(obdata, newdata, all = TRUE))
    #   # p_reform_bar2_terms_full <- terms(rbind(p_reform_bar2_mf, p_reform_bar2_mf_new))
    #   p_reform_bar2_xlev_full <- .getXlevels(p_reform_bar2_terms, rbind(p_reform_bar2_mf, model.frame(reform_bar2, newdata)))
    #   if (!identical(p_reform_bar2_xlev, p_reform_bar2_xlev_full)) {
    #     p_reform_bar2_xlev <- p_reform_bar2_xlev_full
    #   }
    #
    #   partition_index_obdata <- list(reform_bar2_vals = p_reform_bar2_vals, reform_bar2_xlev = p_reform_bar2_xlev)
    #   # partition_index_obdata <- as.vector(model.matrix(reform_bar2, obdata))
    # } else {
    #   reform_bar2 <- NULL
    #   partition_index_obdata <- NULL
    # }

    # matrix cholesky
    if (local_list$method == "all") {
      # if (!is.null(object$random)) {
      #   randcov_names <- get_randcov_names(object$random)
      #   randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
      # }
      # partition_matrix_val <- partition_matrix(object$partition_factor, obdata)
      # cov_matrix_val <- cov_matrix(
      #   spcov_params_val, spdist(obdata, xcoord, ycoord), randcov_params_val,
      #   randcov_Zs, partition_matrix_val
      # )
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

    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      pred_splm <- parallel::parLapply(cl, newdata_list, get_pred_splm,
        se.fit = se.fit,
        interval = interval, formula = object$terms,
        obdata = obdata, xcoord = xcoord, ycoord = ycoord,
        spcov_params_val = spcov_params_val, random = object$random,
        randcov_params_val = randcov_params_val,
        reform_bar2_list = reform_bar2_list,
        Z_index_obdata_list = Z_index_obdata_list,
        reform_bar1_list = reform_bar1_list,
        Z_val_obdata_list = Z_val_obdata_list,
        partition_factor = object$partition_factor,
        reform_bar2 = reform_bar2, partition_index_obdata = partition_index_obdata,
        cov_lowchol = cov_lowchol,
        Xmat = model.matrix(object),
        y = model.response(model.frame(object)),
        offset = model.offset(model.frame(object)), dim_coords = object$dim_coords,
        betahat = coefficients(object), cov_betahat = vcov(object),
        contrasts = object$contrasts,
        local = local_list, xlevels = object$xlevels, diagtol = object$diagtol
      )
      cl <- parallel::stopCluster(cl)
    } else {
      pred_splm <- lapply(newdata_list, get_pred_splm,
        se.fit = se.fit,
        interval = interval, formula = object$terms,
        obdata = obdata, xcoord = xcoord, ycoord = ycoord,
        spcov_params_val = spcov_params_val, random = object$random,
        randcov_params_val = randcov_params_val,
        reform_bar2_list = reform_bar2_list,
        Z_index_obdata_list = Z_index_obdata_list,
        reform_bar1_list = reform_bar1_list,
        Z_val_obdata_list = Z_val_obdata_list,
        partition_factor = object$partition_factor,
        reform_bar2 = reform_bar2, partition_index_obdata = partition_index_obdata,
        cov_lowchol = cov_lowchol,
        Xmat = model.matrix(object),
        y = model.response(model.frame(object)),
        offset = model.offset(model.frame(object)), dim_coords = object$dim_coords,
        betahat = coefficients(object), cov_betahat = vcov(object),
        contrasts = object$contrasts,
        local = local_list, xlevels = object$xlevels, diagtol = object$diagtol
      )
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
        if (add_newdata_rows) {
          names(fit) <- object$missing_index
          names(se) <- object$missing_index
        }
        return(list(fit = fit, se.fit = se))
      } else {
        if (add_newdata_rows) {
          names(fit) <- object$missing_index
        }
        return(fit)
      }
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
      # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
      # tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      fit <- cbind(fit, lwr, upr)
      row.names(fit) <- 1:NROW(fit)
      if (se.fit) {
        if (add_newdata_rows) {
          row.names(fit) <- object$missing_index
          names(se) <- object$missing_index
        }
        return(list(fit = fit, se.fit = se))
      } else {
        if (add_newdata_rows) {
          row.names(fit) <- object$missing_index
        }
        return(fit)
      }
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
    fit <- cbind(fit, lwr, upr)
    row.names(fit) <- 1:NROW(fit)
    if (se.fit) {
      if (add_newdata_rows) {
        row.names(fit) <- object$missing_index
        names(se) <- object$missing_index
      }
      return(list(fit = fit, se.fit = se))
    } else {
      if (add_newdata_rows) {
        row.names(fit) <- object$missing_index
      }
      return(fit)
    }
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}

#' @rdname predict.spmodel
#' @method predict spautor
#' @order 2
#' @export
predict.spautor <- function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf, interval = c("none", "confidence", "prediction"),
                            level = 0.95, type = c("response", "terms"), local, terms = NULL, ...) {

  # match interval argument so the three display
  interval <- match.arg(interval)
  type <- match.arg(type)

  # deal with local
  if (missing(local)) {
    local <- NULL
  }

  # deal with local
  if (is.null(local)) {
    local <- FALSE
  }

  # check scale is numeric (if specified)
  if (!is.null(scale) && !is.numeric(scale)) {
    stop("scale must be numeric.", call. = FALSE)
  }

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }

  # write newdata if predicting missing data
  newdata <- object$data[object$missing_index, , drop = FALSE]

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save randcov param vector
  randcov_params_val <- coef(object, type = "randcov")



  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  if (any(grepl("nmatrix.", attributes(formula_newdata)$dataClasses, fixed = TRUE)) && NROW(newdata) == 1) {
    newdata <- newdata[c(1, 1), , drop = FALSE]
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
    newdata_model <- newdata_model[1, , drop = FALSE]
    # find offset
    offset <- model.offset(newdata_model_frame)
    if (!is.null(offset)) {
      offset <- offset[1]
    }
    newdata <- newdata[1, , drop = FALSE]
  } else {
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    # assumes that predicted observations are not outside the factor levels
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
    # find offset
    offset <- model.offset(newdata_model_frame)
  }
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

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

    # # randcov
    randcov_Zs_val <- get_randcov_Zs(randcov_names = names(randcov_params_val), data = object$data)
    # making the partition matrix
    partition_matrix_val <- partition_matrix(object$partition_factor, object$data)
    # making the covariance matrix
    cov_matrix_val <- cov_matrix(spcov_params_val, object$W, randcov_params_val, randcov_Zs_val, partition_matrix_val, object$M)
    # cov_matrix_val_obs <- covmatrix(object)

    # making the covariance vector
    cov_vector_val <- cov_matrix_val[object$missing_index, object$observed_index, drop = FALSE]
    # cov_vector_val <- covmatrix(object, newdata = object$newdata)

    # splitting the covariance vector
    cov_vector_val_list <- split(cov_vector_val, seq_len(NROW(cov_vector_val)))

    # lower triangular cholesky
    cov_matrix_lowchol <- t(chol(cov_matrix_val[object$observed_index, object$observed_index, drop = FALSE]))
    # cov_matrix_lowchol <- t(chol(cov_matrix_val_obs))

    # find X observed
    X <- model.matrix(object)
    SqrtSigInv_X <- forwardsolve(cov_matrix_lowchol, X)

    # beta hat
    betahat <- coef(object)
    # residuals pearson
    residuals_pearson <- residuals(object, type = "pearson")
    # cov beta hat
    cov_betahat <- vcov(object)
    # total var
    total_var_list <- as.list(diag(cov_matrix_val[object$missing_index, object$missing_index, drop = FALSE]))


    # local prediction list (only for parallel)
    local_list <- get_local_list_prediction(local)

    # local stuff for parallel
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cluster_list <- lapply(seq_along(newdata_model_list), function(l) {
        cluster_list_element <- list(
          x0 = newdata_model_list[[l]],
          c0 = cov_vector_val_list[[l]],
          s0 = total_var_list[[l]]
        )
      })
      pred_spautor <- parallel::parLapply(cl, cluster_list, get_pred_spautor_parallel,
        cov_matrix_lowchol, betahat,
        residuals_pearson,
        cov_betahat, SqrtSigInv_X,
        se.fit = se.fit,
        interval = interval
      )
      cl <- parallel::stopCluster(cl)
    } else {
      # make predictions
      pred_spautor <- mapply(
        x0 = newdata_model_list, c0 = cov_vector_val_list, s0 = total_var_list,
        FUN = function(x0, c0, s0) {
          get_pred_spautor(
            x0 = x0, c0 = c0, s0 = s0,
            cov_matrix_lowchol, betahat,
            residuals_pearson,
            cov_betahat, SqrtSigInv_X,
            se.fit = se.fit,
            interval = interval
          )
        }, SIMPLIFY = FALSE
      )
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
        names(fit) <- object$missing_index
        names(se) <- object$missing_index
        return(list(fit = fit, se.fit = se))
      } else {
        names(fit) <- object$missing_index
        return(fit)
      }
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
      fit <- cbind(fit, lwr, upr)
      row.names(fit) <- 1:NROW(fit)
      if (se.fit) {
        row.names(fit) <- object$missing_index
        names(se) <- object$missing_index
        return(list(fit = fit, se.fit = se))
      } else {
        row.names(fit) <- object$missing_index
        return(fit)
      }
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
    fit <- cbind(fit, lwr, upr)
    row.names(fit) <- 1:NROW(fit)
    if (se.fit) {
      row.names(fit) <- object$missing_index
      names(se) <- object$missing_index
      return(list(fit = fit, se.fit = se))
    } else {
      row.names(fit) <- object$missing_index
      return(fit)
    }
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}

get_pred_splm <- function(newdata_list, se.fit, interval, formula, obdata, xcoord, ycoord,
                          spcov_params_val, random, randcov_params_val, reform_bar2_list,
                          Z_index_obdata_list, reform_bar1_list, Z_val_obdata_list, partition_factor,
                          reform_bar2, partition_index_obdata, cov_lowchol,
                          Xmat, y, offset, betahat, cov_betahat, dim_coords, contrasts, local, xlevels, diagtol = diagtol) {




  # storing partition vector
  partition_vector <- partition_vector(partition_factor,
    data = obdata,
    newdata = newdata_list$row, reform_bar2 = reform_bar2,
    partition_index_data = partition_index_obdata
  )

  # subsetting partition vector (efficient but causes problems later with
  # random effect subsetting)
  if (!is.null(partition_vector) && local$method %in% c("distance", "covariance") &&
    (is.null(random) || !labels(terms(partition_factor)) %in% labels(terms(random)))) {
    partition_index <- as.vector(partition_vector) == 1
    Z_index_obdata_list <- lapply(Z_index_obdata_list, function(x) {
      x$reform_bar2_vals <- x$reform_bar2_vals[partition_index]
      x
    })
    obdata <- obdata[partition_index, , drop = FALSE]
    partition_vector <- Matrix(1, nrow = 1, ncol = NROW(obdata))
  }

  dist_vector <- spdist_vectors(newdata_list$row, obdata, xcoord, ycoord, dim_coords)

  # # subsetting data if method distance
  # if (local$method == "distance") {
  #   n <- length(dist_vector)
  #   nn_index <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, local$size))]
  #   obdata <- obdata[nn_index, , drop = FALSE]
  #   dist_vector <- dist_vector[, nn_index]
  # }

  # making random vector if necessary
  if (!is.null(randcov_params_val)) {
    randcov_vector_val <- randcov_vector(randcov_params_val, obdata, newdata_list$row, reform_bar2_list, Z_index_obdata_list)
  } else {
    randcov_vector_val <- NULL
  }

  # making the covariance vector
  cov_vector_val <- cov_vector(spcov_params_val, dist_vector, randcov_vector_val, partition_vector)

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
    cov_index <- order(as.numeric(cov_vector_val))[seq(from = n, to = max(1, n - local$size + 1))] # use abs() here?
    obdata <- obdata[cov_index, , drop = FALSE]
    cov_vector_val <- cov_vector_val[cov_index]
  }

  if (local$method %in% c("distance", "covariance")) {
    if (!is.null(random)) {
      randcov_names <- get_randcov_names(random)
      xlev_list <- lapply(Z_index_obdata_list, function(x) x$reform_bar2_xlev)
      randcov_Zs <- get_randcov_Zs(obdata, randcov_names, xlev_list = xlev_list)
    }
    partition_matrix_val <- partition_matrix(partition_factor, obdata)
    cov_matrix_val <- cov_matrix(
      spcov_params_val, spdist(obdata, xcoord, ycoord), randcov_params_val,
      randcov_Zs, partition_matrix_val,
      diagtol = diagtol
    )
    cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(cov_matrix_val)))
    model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass, xlev = xlevels)
    Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
    y <- model.response(model_frame)
    offset <- model.offset(model_frame)
  }


  # handle offset
  if (!is.null(offset)) {
    y <- y - offset
  }

  c0 <- as.numeric(cov_vector_val)
  SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
  SqrtSigInv_y <- forwardsolve(cov_lowchol, y)
  residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
  x0 <- newdata_list$x0

  fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
  H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  if (se.fit || interval == "prediction") {
    total_var <- sum(spcov_params_val[["de"]], spcov_params_val[["ie"]], randcov_params_val)
    var <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
    pred_list <- list(fit = fit, var = var)
  } else {
    pred_list <- list(fit = fit)
  }
  pred_list
}

get_pred_spautor <- function(x0, c0, s0, cov_matrix_lowchol, betahat, residuals_pearson, cov_betahat, SqrtSigInv_X, se.fit, interval) {
  SqrtSigInv_c0 <- forwardsolve(cov_matrix_lowchol, c0)
  fit <- as.numeric(x0 %*% betahat + crossprod(SqrtSigInv_c0, residuals_pearson))
  if (se.fit || interval == "prediction") {
    H <- x0 - crossprod(SqrtSigInv_c0, SqrtSigInv_X)
    var <- as.numeric(s0 - crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% tcrossprod(cov_betahat, H))
    pred_list <- list(fit = fit, var = var)
  } else {
    pred_list <- list(fit = fit)
  }
  pred_list
}

get_pred_spautor_parallel <- function(cluster_list, cov_matrix_lowchol, betahat, residuals_pearson, cov_betahat, SqrtSigInv_X, se.fit, interval) {
  x0 <- cluster_list$x0
  c0 <- cluster_list$c0
  s0 <- cluster_list$s0
  get_pred_spautor(x0, c0, s0, cov_matrix_lowchol, betahat, residuals_pearson, cov_betahat, SqrtSigInv_X, se.fit, interval)
}


#' @name predict.spmodel
#' @method predict splm_list
#' @order 3
#' @export
predict.splm_list <- function(object, newdata, se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                              level = 0.95, local, ...) {
  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with local
  if (missing(local)) {
    local <- NULL
  }

  if (missing(newdata)) {
    preds <- lapply(object, function(x) {
      predict(x, se.fit = se.fit, interval = interval, level = level, local = local, ...)
    })
  } else {
    preds <- lapply(object, function(x) {
      predict(x, newdata = newdata, se.fit = se.fit, interval = interval, level = level, local = local, ...)
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
    if (missing(local)) {
      local <- NULL
    }
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
    if (missing(local)) {
      local <- NULL
    }
    # do spautor prediction
    spautor_pred <- do.call(predict, list(object = object$spautor, newdata = newdata, local = local))
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
    if (missing(local)) {
      local <- NULL
    }
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
    if (missing(local)) {
      local <- NULL
    }
    # do spautor prediction
    spautor_pred <- lapply(object$spautor_list, function(x) {
      do.call(predict, list(object = x, newdata = newdata, local = local))
    })
  }
  # obtain final predictions
  sprf_pred <- lapply(spautor_pred, function(x) ranger_pred$predictions + x)
  names(sprf_pred) <- names(object$spautor_list)
  sprf_pred
}


# extra lists for use with loocv (do nothing but clean with predict)
get_extra_randcov_list <- function(object, obdata, newdata) {
  # random stuff
  if (!is.null(object$random)) {
    randcov_names <- get_randcov_names(object$random)
    # this causes a memory leak and was not even needed
    # randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
    # comment out here for simple
    reform_bar_list <- lapply(randcov_names, function(randcov_name) {
      bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
      reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
      if (bar_split[[1]] != "1") {
        reform_bar1 <- reformulate(bar_split[[1]], intercept = FALSE)
      } else {
        reform_bar1 <- NULL
      }
      list(reform_bar2 = reform_bar2, reform_bar1 = reform_bar1)
    })
    reform_bar2_list <- lapply(reform_bar_list, function(x) x$reform_bar2)
    names(reform_bar2_list) <- randcov_names
    reform_bar1_list <- lapply(reform_bar_list, function(x) x$reform_bar1)
    names(reform_bar1_list) <- randcov_names
    Z_index_obdata_list <- lapply(reform_bar2_list, function(reform_bar2) {

      reform_bar2_mf <- model.frame(reform_bar2, obdata)
      reform_bar2_terms <- terms(reform_bar2_mf)
      reform_bar2_xlev <- .getXlevels(reform_bar2_terms, reform_bar2_mf)
      reform_bar2_mx <- model.matrix(reform_bar2, obdata)
      reform_bar2_names <- colnames(reform_bar2_mx)
      reform_bar2_split <- split(reform_bar2_mx, seq_len(NROW(reform_bar2_mx)))
      reform_bar2_vals <- reform_bar2_names[vapply(reform_bar2_split, function(y) which(as.logical(y)), numeric(1))]

      # adding dummy levels if newdata observations of random effects are not in original data
      # terms object is unchanged if levels change
      # reform_bar2_mf_new <- model.frame(reform_bar2, newdata)
      # reform_bar2_mf_full <- model.frame(reform_bar2, merge(obdata, newdata, all = TRUE))
      # reform_bar2_terms_full <- terms(rbind(reform_bar2_mf, reform_bar2_mf_new))
      reform_bar2_xlev_full <- .getXlevels(reform_bar2_terms, rbind(reform_bar2_mf, model.frame(reform_bar2, newdata)))
      if (!identical(reform_bar2_xlev, reform_bar2_xlev_full)) {
        reform_bar2_xlev <- reform_bar2_xlev_full
      }

      list(reform_bar2_vals = reform_bar2_vals, reform_bar2_xlev = reform_bar2_xlev)
    })
    # Z_index_obdata_list <- lapply(reform_bar2_list, function(reform_bar2) as.vector(model.matrix(reform_bar2, obdata)))
    names(Z_index_obdata_list) <- randcov_names
    Z_val_obdata_list <- lapply(reform_bar1_list, function(reform_bar1) {
      if (is.null(reform_bar1)) {
        return(NULL)
      } else {
        return(as.vector(model.matrix(reform_bar1, obdata)))
      }
    })
    names(Z_val_obdata_list) <- randcov_names
  } else {
    reform_bar2_list <- NULL
    Z_index_obdata_list <- NULL
    reform_bar1_list <- NULL
    Z_val_obdata_list <- NULL
  }
  list(reform_bar1_list = reform_bar1_list, reform_bar2_list = reform_bar2_list,
       Z_val_obdata_list = Z_val_obdata_list, Z_index_obdata_list = Z_index_obdata_list)
}

get_extra_partition_list <- function(object, obdata, newdata) {
  # partition factor stuff
  if (!is.null(object$partition_factor)) {
    partition_factor_val <- get_partition_name(labels(terms(object$partition_factor)))
    bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
    p_reform_bar2_mf <- model.frame(reform_bar2, obdata)
    p_reform_bar2_terms <- terms(p_reform_bar2_mf)
    p_reform_bar2_xlev <- .getXlevels(p_reform_bar2_terms, p_reform_bar2_mf)
    p_reform_bar2_mx <- model.matrix(reform_bar2, obdata)
    p_reform_bar2_names <- colnames(p_reform_bar2_mx)
    p_reform_bar2_split <- split(p_reform_bar2_mx, seq_len(NROW(p_reform_bar2_mx)))
    p_reform_bar2_vals <- p_reform_bar2_names[vapply(p_reform_bar2_split, function(y) which(as.logical(y)), numeric(1))]


    # adding dummy levels if newdata observations of random effects are not in original data
    # terms object is unchanged if levels change
    # p_reform_bar2_mf_new <- model.frame(reform_bar2, newdata)
    # reform_bar2_mf_full <- model.frame(reform_bar2, merge(obdata, newdata, all = TRUE))
    # p_reform_bar2_terms_full <- terms(rbind(p_reform_bar2_mf, p_reform_bar2_mf_new))
    p_reform_bar2_xlev_full <- .getXlevels(p_reform_bar2_terms, rbind(p_reform_bar2_mf, model.frame(reform_bar2, newdata)))
    if (!identical(p_reform_bar2_xlev, p_reform_bar2_xlev_full)) {
      p_reform_bar2_xlev <- p_reform_bar2_xlev_full
    }

    partition_index_obdata <- list(reform_bar2_vals = p_reform_bar2_vals, reform_bar2_xlev = p_reform_bar2_xlev)
    # partition_index_obdata <- as.vector(model.matrix(reform_bar2, obdata))
  } else {
    reform_bar2 <- NULL
    partition_index_obdata <- NULL
  }
  list(reform_bar2 = reform_bar2, partition_index_obdata = partition_index_obdata)
}
