#' Model predictions (Kriging)
#'
#' @description Predicted values and intervals based on a fitted model object.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param newdata A data frame or \code{sf} object in which to
#'   look for variables with which to predict. If a data frame, \code{newdata}
#'   must contain all variables used by \code{formula(object)} and all variables
#'   representing coordinates. If an \code{sf} object, \code{newdata} must contain
#'   all variables used by \code{formula(object)} and coordinates are obtained
#'   from the geometry of \code{newdata}. If omitted, missing data from the
#'   fitted model object are used.
#' @param se.fit A logical indicating if standard errors are returned.
#'   The default is \code{FALSE}.
#' @param interval Type of interval calculation. The default is \code{"none"}.
#'   Other options are \code{"confidence"} (for confidence intervals) and
#'   \code{"prediction"} (for prediction intervals).
#' @param level Tolerance/confidence level. The default is \code{0.95}.
#' @param local A optional logical or list controlling the big data approximation. If omitted, \code{local}
#'   is set to \code{TRUE} or \code{FALSE} based on the sample size of the fitted
#'   model object and/or the prediction size of \code{newdata} -- if the sample
#'   size or prediction size exceeds 5000, \code{local} is
#'   set to \code{TRUE}, otherwise it is set to \code{FALSE}. If \code{FALSE}, no big data approximation
#'   is implemented. If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item{\code{method}: }{The big data approximation method. If \code{method = "all"},
#'       all observations are used and \code{size} is ignored. If \code{method = "distance"},
#'       the \code{size} data observations closest (in terms of Euclidean distance)
#'       to the observation requiring prediction are used.
#'       If \code{method = "covariance"}, the \code{size} data observations
#'       with the highest covariance with the observation requiring prediction are used.
#'       If random effects and partition factors are not used in estimation and
#'       the spatial covariance function is monotone decreasing,
#'       \code{"distance"} and \code{"covariance"} are equivalent. The default
#'       is \code{"covariance"}. Only used with models fit using [splm()].}
#'     \item{\code{size}: }{The number of data observations to use when \code{method}
#'       is \code{"distance"} or \code{"covariance"}. The default is 50. Only used
#'       with models fit using [splm()].}
#'     \item{\code{parallel}: }{If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. The default is \code{FALSE}.}
#'     \item{\code{ncores}: }{If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.}
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 50, method = "covariance", parallel = FALSE)}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details The (empirical) best linear unbiased predictions (i.e., Kriging
#'   predictions) at each site are returned when \code{interval} is \code{"none"}
#'   or \code{"prediction"} alongside standard errors. Prediction intervals
#'   are also returned if \code{interval} is \code{"prediction"}. When
#'   \code{interval} is \code{"confidence"}, the estimated mean is returned
#'   alongside standard errors and confidence intervals for the mean.
#'
#' @return If \code{se.fit} is \code{FALSE}, \code{predict.spmod()} returns
#'   a vector of predictions or a matrix of predictions with column names
#'   \code{fit}, \code{lwr}, and \code{upr} if \code{interval} is \code{"confidence"}
#'   or \code{"prediction"}.
#'
#'   If \code{se.fit} is \code{TRUE}, a list with the following components is returned:
#'   \itemize{
#'     \item{\code{fit}: }{vector or matrix as above}
#'     \item{\code{se.fit: }}{standard error of each fit}
#'   }
#'
#' @method predict spmod
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
predict.spmod <- function(object, newdata, se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                          level = 0.95, local, ...) {

  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with local
  if (missing(local)) {
    local <- NULL
  }

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("newdata must be specified in the newdata argument or in the newdata element of object.", call. = FALSE)
  }

  switch(object$fn,
    "splm" = predict_splm(object, newdata, se.fit, interval, level, local, ...),
    "spautor" = predict_spautor(object, se.fit, interval, level, local, ...)
  )
}
predict_splm <- function(object, newdata, se.fit = FALSE, interval,
                         level = 0.95, local, ...) {

  # browser()

  # rename relevant quantities
  obdata <- object$obdata
  xcoord <- object$xcoord
  ycoord <- object$ycoord

  # write newdata if predicting missing data
  if (missing(newdata)) {
    add_newdata_rows <- TRUE
    if (is.null(object$newdata)) {
      # stop("No missing data to predict", call. = FALSE)
    }
    newdata <- object$newdata
  } else {
    add_newdata_rows <- FALSE
  }

  # deal with local
  if (is.null(local)) {
    if (object$n > 5000 || NROW(newdata) > 5000) {
      local <- TRUE
      message("Because either the sample size of the fitted model object or the number of desired predictions exceeds 5000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun predict() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save randcov param vector
  randcov_params_val <- coef(object, type = "randcov")

  # if (inherits(newdata, "SpatialPointsDataFrame")) { # may want to handle this better and ensure data was sf object
  #   newdata <- sp_to_df(newdata)
  #   names(newdata)[[which(names(newdata) == "xcoord")]] <- xcoord # only relevant if newdata is sp data is not
  #   names(newdata)[[which(names(newdata) == "ycoord")]] <- ycoord # only relevant if newdata is sp data is not
  # }

  attr_sp <- attr(class(newdata), "package")
  if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
    # if (inherits(newdata, c("SpatialPointsDataFrame", "SpatialPolygonsDataFrame"))) {
    # if (!requireNamespace("sf", quietly = TRUE)) { # requireNamespace checks if sf is installed
    #   stop("Install the sf R package to use sp objects in splm()", call. = FALSE)
    # } else {
    #   newdata <- sf::st_as_sf(newdata)
    # }
    # newdata <- sf::st_as_sf(newdata)
    stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
  }

  if (inherits(newdata, "sf")) {

    # turn polygon into centroid
    # if (!requireNamespace("sf", quietly = TRUE)) { # requireNamespace checks if sf is installed
    #   stop("Install the sf R package to use the centroid of sf POLYGON objects in splm()", call. = FALSE)
    # } else {
    #   # warning here from sf about assuming attributes are constant over geometries of x
    #   newdata <- suppressWarnings(sf::st_centroid(newdata))
    # }
    newdata <- suppressWarnings(sf::st_centroid(newdata))

    newdata <- sf_to_df(newdata)
    names(newdata)[[which(names(newdata) == "xcoord")]] <- as.character(xcoord) # only relevant if newdata is sf data is not
    names(newdata)[[which(names(newdata) == "ycoord")]] <- as.character(ycoord) # only relevant if newdata is sf data is not
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

  # storing new model matrix
  # labels_newdata <- labels(terms(object$formula))
  # if (length(labels_newdata) == 0) {
  #   labels_newdata <- "1"
  # }
  # formula_newdata <- reformulate(labels_newdata)
  formula_newdata <- update(formula(object), NULL ~ .)
  newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass)
  # assumes that predicted observations are not outside the factor levels
  newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # storing newdata as a list
  newdata_list <- split(newdata, seq_len(NROW(newdata)))



  if (interval %in% c("none", "prediction")) {

    # local prediction list
    local_list <- get_local_list_prediction(local)

    # hardcoding none covariance
    if (inherits(spcov_params_val, "none") && is.null(randcov_params_val)) {
      local_list$method <- "all"
    }

    # random stuff
    if (!is.null(object$random)) {
      randcov_names <- get_randcov_names(object$random)
      randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
      # comment out here for simple
      reform_bar_list <- lapply(randcov_names, function(randcov_name) {
        bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
        reform_bar2 <- reformulate(paste0("as.numeric(", bar_split[[2]], ")"), intercept = FALSE)
        if (bar_split[[1]] != "1") {
          reform_bar1 <- reformulate(paste0("as.numeric(", bar_split[[1]], ")"), intercept = FALSE)
        } else {
          reform_bar1 <- NULL
        }
        list(reform_bar2 = reform_bar2, reform_bar1 = reform_bar1)
      })
      reform_bar2_list <- lapply(reform_bar_list, function(x) x$reform_bar2)
      names(reform_bar2_list) <- randcov_names
      reform_bar1_list <- lapply(reform_bar_list, function(x) x$reform_bar1)
      names(reform_bar1_list) <- randcov_names
      Z_index_obdata_list <- lapply(reform_bar2_list, function(reform_bar2) as.vector(model.matrix(reform_bar2, obdata)))
      names(Z_index_obdata_list) <- randcov_names
      Z_val_obdata_list <- lapply(reform_bar1_list, function(reform_bar1) {
        if (is.null(reform_bar1)) {
          return(NULL)
        } else {
          return(as.vector(model.matrix(reform_bar1, obdata)))
        }
      })
      names(Z_val_obdata_list) <- randcov_names
      # # uncomment for simple
      # reform_bar2_list <- NULL
      # Z_index_obdata_list <- NULL
      # reform_bar1_list <- NULL
      # Z_val_obdata_list <- NULL
    } else {
      reform_bar2_list <- NULL
      Z_index_obdata_list <- NULL
      reform_bar1_list <- NULL
      Z_val_obdata_list <- NULL
    }

    # partition factor stuff
    if (!is.null(object$partition_factor)) {
      # comment out for simple
      # browser()
      partition_factor_val <- get_partition_name(labels(terms(object$partition_factor)))
      bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
      reform_bar2 <- reformulate(paste0("as.numeric(", bar_split[[2]], ")"), intercept = FALSE)
      partition_index_obdata <- as.vector(model.matrix(reform_bar2, obdata))
      # # uncomment for simple
      # reform_bar2 <- NULL
      # partition_index_obdata <- NULL
    } else {
      reform_bar2 <- NULL
      partition_index_obdata <- NULL
    }

    # matrix cholesky
    if (local_list$method == "all") {
      if (!is.null(object$random)) {
        randcov_names <- get_randcov_names(object$random)
        randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
      }
      partition_matrix_val <- partition_matrix(object$partition_factor, obdata)
      cov_matrix_val <- cov_matrix(
        spcov_params_val, spdist(obdata, xcoord, ycoord), randcov_params_val,
        randcov_Zs, partition_matrix_val
      )
      cov_lowchol <- t(chol(cov_matrix_val))
    } else {
      cov_lowchol <- NULL
    }

    # browser()
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      pred_splm <- parallel::parLapply(cl, newdata_list, get_pred_splm,
        se.fit = se.fit,
        interval = interval, formula = object$formula,
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
        y = model.response(model.frame(object)), dim_coords = object$dim_coords,
        betahat = coefficients(object), cov_betahat = vcov(object),
        contrasts = object$contrasts,
        local = local_list
      )
      cl <- parallel::stopCluster(cl)
    } else {
      pred_splm <- lapply(newdata_list, get_pred_splm,
        se.fit = se.fit,
        interval = interval, formula = object$formula,
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
        y = model.response(model.frame(object)), dim_coords = object$dim_coords,
        betahat = coefficients(object), cov_betahat = vcov(object),
        contrasts = object$contrasts,
        local = local_list
      )
    }

    if (interval == "none") {
      fit <- vapply(pred_splm, function(x) x$fit, numeric(1))
      if (se.fit) {
        vars <- vapply(pred_splm, function(x) x$var, numeric(1))
        se <- sqrt(vars)
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
      vars <- vapply(pred_splm, function(x) x$var, numeric(1))
      se <- sqrt(vars)
      # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
      tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      fit <- cbind(fit, lwr, upr)
      rownames(fit) <- 1:NROW(fit)
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
    newdata_model_list <- split(newdata_model, seq_len(NROW(newdata_model)))
    vars <- as.numeric(vapply(newdata_model_list, function(x) crossprod(x, vcov(object) %*% x), numeric(1)))
    se <- sqrt(vars)
    # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
    tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    fit <- cbind(fit, lwr, upr)
    # rownames(fit) <- 1:NROW(fit)
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

#
#
#
#     # splitting up new data model matrix into list
#     newdata_model_split <- split(newdata_model, seq_len(NROW(newdata_model)))
#
#     # making distance vector
#     ## finding 1D distance
#     if (object$dim_coords == 1) {
#       dist_vector <- sqrt(outer(X = newdata[[object$xcoord]], Y = object$data[[object$xcoord]], FUN = function(X, Y) (X - Y)^2))
#     } else if (object$dim_coord == 2) { ## finding 2D distance
#       dist_vector_x <- outer(X = newdata[[object$xcoord]], Y = object$data[[object$xcoord]], FUN = function(X, Y) (X - Y)^2)
#       dist_vector_y <- outer(X = newdata[[object$ycoord]], Y = object$data[[object$ycoord]], FUN = function(X, Y) (X - Y)^2)
#       dist_vector <- sqrt(dist_vector_x + dist_vector_y)
#     } else {
#       dist_vector <- matrix(Inf, nrow = NROW(newdata), ncol = NROW(object$data)) ## 0D distance (coords not used for "none")
#     }
#     #browser()
#     # making partition_vector
#     partition_vector <- partition_vector(object$partition_factor, object$data, newdata)
#
#     # making the covariance vector
#     cov_vector_val <- cov_vector(spcov_params_val, dist_vector,
#                                  coef(object, type = "randcov"),
#                                  partition_vector)
#
#     cov_vector_val_split <- split(cov_vector_val, seq_len(NROW(cov_vector_val)))
#
#     ## STOPPED HERE 3/26 @ 8:30 AM
#
#   }
#
#   if (interval == "none") {
#     # just returning fitted values
#     newdata_model_split <- split(newdata_model, 1:NROW(newdata_model))
#     if (object$data_object$dim_coords == 1) {
#       dist_vector <- sqrt(outer(X = newdata_xcoord, Y = data_xcoord, FUN = function(X, Y) (X - Y)^2))
#     } else if (object$data_object$dim_coords == 2) {
#       dist_vector_x <- outer(X = newdata_xcoord, Y = data_xcoord, FUN = function(X, Y) (X - Y)^2)
#       dist_vector_y <- outer(X = newdata_ycoord, Y = data_ycoord, FUN = function(X, Y) (X - Y)^2)
#       dist_vector <- sqrt(dist_vector_x + dist_vector_y)
#     } else {
#       dist_vector <- matrix(Inf, nrow = NROW(newdata), ncol = NROW(object$data))
#     }
#
#
#     # making the covariance vector
#     cov_vector_val <- cov_vector(spcov_params_val, dist_vector,
#                                  coef(object, type = "randcov"),
#                                  object$partition_factor, object$data, newdata)
#
#     cov_vector_val_split <- split(cov_vector_val, seq_len(NROW(cov_vector_val)))
#
#     # making a covariance matrix
#     if (length(object$coords) == 1) {
#       dist_matrix <- spdist(data_xcoord)
#     } else if (length(object$coords) == 2) {
#       dist_matrix <- spdist(data_xcoord, data_ycoord)
#     } else {
#       dist_matrix <- matrix(Inf, nrow = NROW(object$data), NCOL(object$data))
#       diag(dist_matrix) <- 0
#     }
#
#
#
#     # randcov Zs
#     randcov_Zs_val <- get_randcov_Zs(names(randcov_params_val), object$data)
#     # partition matrix
#     partition_matrix_val <- partition_matrix(object$partition_factor, object$data)
#
#     # making the covariance matrix
#     cov_matrix_val <- cov_matrix(spcov_params_val, dist_matrix, randcov_params_val, randcov_Zs_val, partition_matrix_val)
#
#     # lower triangular cholesky
#     cov_matrix_lowchol <- t(chol(cov_matrix_val))
#     SqrtSigInv_X <- forwardsolve(cov_matrix_lowchol, object$model$X)
#
#     # beta hat
#     betahat <- coef(object)
#     # residuals pearson / pearson
#     residuals_pearson <- residuals(object, type = "pearson")
#     # cov beta hat
#     cov_betahat <- vcov(object)
#     # total var
#     total_var <- sum(spcov_params_val[["de"]], spcov_params_val[["ie"]], randcov_params_val)
#     preds_val <- mapply(x0 = newdata_model_split, c0 = cov_vector_val_split,
#                         FUN = function(x0, c0) get_preds_splm(x0 = x0, c0 = c0,
#                                                          cov_matrix_lowchol, betahat,
#                                                          residuals_pearson, total_var,
#                                                          cov_betahat, SqrtSigInv_X), SIMPLIFY = FALSE)
#     fit <- vapply(preds_val, function(x) x$fit, numeric(1))
#     vars <- vapply(preds_val, function(x) x$var, numeric(1))
#     se <- sqrt(vars)
#     if (se.fit) {
#       return(list(fit = fit, se.fit = se))
#     } else {
#       return(fit)
#     }
#   } else if (interval == "confidence") {
#     # finding fitted values of the mean parameters
#     fit <- as.numeric(newdata_model %*% coef(object))
#     newdata_model_split <- split(newdata_model, 1:NROW(newdata_model))
#     vars <- as.numeric(vapply(newdata_model_split, function(x) crossprod(x, vcov(object) %*% x), numeric(1)))
#     se <- sqrt(vars)
#     # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
#     tstar <- qnorm(1 - (1 - level) / 2)
#     lwr <- fit - tstar * se
#     upr <- fit + tstar * se
#     fit <- cbind(fit, lwr, upr)
#     rownames(fit) <- 1:NROW(fit)
#     if (se.fit) {
#       return(list(fit = fit, se.fit = se))
#     } else {
#       return(fit)
#     }
#   } else if (interval == "prediction") {
#     newdata_model_split <- split(newdata_model, 1:NROW(newdata_model))
#     if (length(object$coords) == 1) {
#       dist_vector <- sqrt(outer(X = newdata_xcoord, Y = data_xcoord, FUN = function(X, Y) (X - Y)^2))
#     } else if (length(object$coords) == 2) {
#       dist_vector_x <- outer(X = newdata_xcoord, Y = data_xcoord, FUN = function(X, Y) (X - Y)^2)
#       dist_vector_y <- outer(X = newdata_ycoord, Y = data_ycoord, FUN = function(X, Y) (X - Y)^2)
#       dist_vector <- sqrt(dist_vector_x + dist_vector_y)
#     } else {
#       dist_vector <- matrix(Inf, nrow = NROW(newdata), ncol = NROW(object$data))
#     }
#
#     # making the covariance vector
#     cov_vector_val <- cov_vector(spcov_params_val, dist_vector,
#                                  coef(object, type = "randcov"),
#                                  object$partition_factor, object$data, newdata)
#
#     cov_vector_val_split <- split(cov_vector_val, seq_len(NROW(cov_vector_val)))
#
#     # making a covariance matrix
#     if (length(object$coords) == 1) {
#       dist_matrix <- spdist(data_xcoord)
#     } else if (length(object$coords) == 2) {
#       dist_matrix <- spdist(data_xcoord, data_ycoord)
#     } else {
#       dist_matrix <- matrix(Inf, nrow = NROW(object$data), NCOL(object$data))
#       diag(dist_matrix) <- 0
#     }
#
#     # randcov Zs
#     randcov_Zs_val <- get_randcov_Zs(names(randcov_params_val), object$data)
#     # partition matrix
#     partition_matrix_val <- partition_matrix(object$partition_factor, object$data)
#
#     # making the covariance matrix
#     cov_matrix_val <- cov_matrix(spcov_params_val, dist_matrix, randcov_params_val, randcov_Zs_val, partition_matrix_val)
#
#     # lower triangular cholesky
#     cov_matrix_lowchol <- t(chol(cov_matrix_val))
#     SqrtSigInv_X <- forwardsolve(cov_matrix_lowchol, object$model$X)
#
#     # beta hat
#     betahat <- coef(object)
#     # residuals pearson
#     residuals_pearson <- residuals(object, type = "pearson")
#     # cov beta hat
#     cov_betahat <- vcov(object)
#     # total var
#     total_var <- sum(spcov_params_val[["de"]], spcov_params_val[["ie"]], randcov_params_val)
#     preds_val <- mapply(x0 = newdata_model_split, c0 = cov_vector_val_split,
#                         FUN = function(x0, c0) get_preds_splm(x0 = x0, c0 = c0,
#                                                          cov_matrix_lowchol, betahat,
#                                                          residuals_pearson, total_var,
#                                                          cov_betahat, SqrtSigInv_X), SIMPLIFY = FALSE)
#
#     fit <- vapply(preds_val, function(x) x$fit, numeric(1))
#     vars <- vapply(preds_val, function(x) x$var, numeric(1))
#     se <- sqrt(vars)
#     # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
#     tstar <- qnorm(1 - (1 - level) / 2)
#     lwr <- fit - tstar * se
#     upr <- fit + tstar * se
#     fit <- cbind(fit, lwr, upr)
#     rownames(fit) <- 1:NROW(fit)
#     if (se.fit) {
#       return(list(fit = fit, se.fit = se))
#     } else {
#       return(fit)
#     }
#   } else {
#     stop("interval must be none, confidence, or prediction")
#   }
# }

predict_spautor <- function(object, se.fit = FALSE, interval,
                            level = 0.95, local, ...) {

  # browser()

  # deal with local
  if (is.null(local)) {
    local <- FALSE
  }
  # data obs
  # obdata <- object$obdata

  # write newdata if predicting missing data
  newdata <- object$data[object$missing_index, , drop = FALSE]
  if (NROW(newdata) == 0) {
    stop("No missing data to predict", call. = FALSE)
  }

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save randcov param vector
  randcov_params_val <- coef(object, type = "randcov")



  # storing new model matrix
  # labels_newdata <- labels(terms(object$formula))
  # if (length(labels_newdata) == 0) {
  #   labels_newdata <- "1"
  # }
  # formula_newdata <- reformulate(labels_newdata)
  formula_newdata <- update(formula(object), NULL ~ .)
  newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass)
  # note that this will return an error if value of factor used to predict that was not used to fit
  newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # storing newdata as a list
  newdata_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  if (interval %in% c("none", "prediction")) {

    # randcov
    randcov_Zs_val <- get_randcov_Zs(randcov_names = names(randcov_params_val), data = object$data)
    # making the partition matrix
    partition_matrix_val <- partition_matrix(object$partition_factor, object$data)
    # making the covariance matrix
    cov_matrix_val <- cov_matrix(spcov_params_val, object$W, randcov_params_val, randcov_Zs_val, partition_matrix_val, object$M)

    # making the covariance vector
    cov_vector_val <- cov_matrix_val[object$missing_index, object$observed_index, drop = FALSE]

    # splitting the covariance vector
    cov_vector_val_list <- split(cov_vector_val, seq_len(NROW(cov_vector_val)))

    # lower triangular cholesky
    cov_matrix_lowchol <- t(chol(cov_matrix_val[object$observed_index, object$observed_index, drop = FALSE]))

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
    # browser()
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
      if (se.fit) {
        vars <- vapply(pred_spautor, function(x) x$var, numeric(1))
        se <- sqrt(vars)
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
      vars <- vapply(pred_spautor, function(x) x$var, numeric(1))
      se <- sqrt(vars)
      # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
      tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      fit <- cbind(fit, lwr, upr)
      rownames(fit) <- 1:NROW(fit)
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
    vars <- as.numeric(vapply(newdata_model_list, function(x) crossprod(x, vcov(object) %*% x), numeric(1)))
    se <- sqrt(vars)
    # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
    tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    fit <- cbind(fit, lwr, upr)
    # rownames(fit) <- 1:NROW(fit)
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

# get_preds_splm <- function(x0, c0, cov_matrix_lowchol, betahat, residuals_pearson, total_var, cov_betahat, SqrtSigInv_X) {
#   SqrtSigInv_c0 <- forwardsolve(cov_matrix_lowchol, c0)
#   fit <- as.numeric(x0 %*% betahat + crossprod(SqrtSigInv_c0, residuals_pearson))
#   H <- x0 - crossprod(SqrtSigInv_c0, SqrtSigInv_X)
#   var <- as.numeric(total_var - crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% tcrossprod(cov_betahat, H))
#   list(fit = fit, var = var)
# }
#
# get_preds_spautor <- function(x0, c0, s0, cov_matrix_lowchol, betahat, residuals_pearson, cov_betahat, SqrtSigInv_X) {
#   SqrtSigInv_c0 <- forwardsolve(cov_matrix_lowchol, c0)
#   fit <- as.numeric(x0 %*% betahat + crossprod(SqrtSigInv_c0, residuals_pearson))
#   H <- x0 - crossprod(SqrtSigInv_c0, SqrtSigInv_X)
#   var <- as.numeric(s0 - crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% tcrossprod(cov_betahat, H))
#   list(fit = fit, var = var)
# }





get_pred_splm <- function(newdata_row, se.fit, interval, formula, obdata, xcoord, ycoord,
                          spcov_params_val, random, randcov_params_val, reform_bar2_list,
                          Z_index_obdata_list, reform_bar1_list, Z_val_obdata_list, partition_factor,
                          reform_bar2, partition_index_obdata, cov_lowchol,
                          Xmat, y, betahat, cov_betahat, dim_coords, contrasts, local) {

  # browser()

  # storing partition vector
  partition_vector <- partition_vector(partition_factor,
    data = obdata,
    newdata = newdata_row, reform_bar2 = reform_bar2,
    partition_index_data = partition_index_obdata
  )

  # subsetting partition vector
  if (!is.null(partition_vector) && local$method %in% c("distance", "covariance") &&
    !labels(terms(partition_factor)) %in% labels(terms(random))) {
    obdata <- obdata[as.vector(partition_vector) == 1, , drop = FALSE]
    partition_vector <- Matrix(1, nrow = 1, ncol = NROW(obdata))
  }

  dist_vector <- spdist_vectors(newdata_row, obdata, xcoord, ycoord, dim_coords)

  # subsetting data if method distance
  if (local$method == "distance") {
    # nn_list <- nabor::knn(cbind(obdata[[xcoord]], obdata[[ycoord]]),
    #                             cbind(newdata_row[[xcoord]], newdata_row[[ycoord]]), k = max(local$size, NROW(data)))
    # dist_vector <- nn_list$nn.dists
    # nn_index <- as.vector(nn_list$nn.idx)
    # obdata <- obdata[nn_index, , drop = FALSE]
    n <- length(dist_vector)
    nn_index <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, local$size))]
    obdata <- obdata[nn_index, , drop = FALSE]
    dist_vector <- dist_vector[, nn_index]
  }

  # making random vector if necessary
  if (!is.null(randcov_params_val)) {
    randcov_vector_val <- randcov_vector(randcov_params_val, obdata, newdata_row, reform_bar2_list, Z_index_obdata_list)
  } else {
    randcov_vector_val <- NULL
  }

  # making the covariance vector
  cov_vector_val <- cov_vector(spcov_params_val, dist_vector, randcov_vector_val, partition_vector)


  if (local$method == "covariance") {
    n <- length(cov_vector_val)
    cov_index <- order(as.numeric(cov_vector_val))[seq(from = n, to = max(1, n - local$size + 1))] # use abs() here?
    obdata <- obdata[cov_index, , drop = FALSE]
    cov_vector_val <- cov_vector_val[cov_index]
  }

  if (local$method %in% c("distance", "covariance")) {
    if (!is.null(random)) {
      randcov_names <- get_randcov_names(random)
      randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
    }
    partition_matrix_val <- partition_matrix(partition_factor, obdata)
    cov_matrix_val <- cov_matrix(
      spcov_params_val, spdist(obdata, xcoord, ycoord), randcov_params_val,
      randcov_Zs, partition_matrix_val
    )
    cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(cov_matrix_val)))
    model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass)
    Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
    y <- model.response(model_frame)
  }



  c0 <- as.numeric(cov_vector_val)
  SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
  SqrtSigInv_y <- forwardsolve(cov_lowchol, y)
  residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
  x0 <- model.matrix(formula(delete.response(terms(formula))), newdata_row, contrasts = contrasts)

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
