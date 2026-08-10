#' Predict a block average for an \code{splm()} model
#'
#' @param object A fitted model object from [splm()]
#' @param newdata The data (points within the block) requiring prediction
#' @param se.fit Whether to compute the prediction standard error
#' @param scale A scale multiplier for the standard error (or \code{NULL})
#' @param df Degrees of freedom used for interval quantiles
#' @param interval The type of interval (\code{"none"}, \code{"confidence"}, or \code{"prediction"})
#' @param level The confidence/prediction interval level
#' @param type The prediction type (only \code{"response"} and \code{"terms"} are relevant here)
#' @param local A list or logical controlling the big data approximation
#' @param terms Terms to use when \code{type} is \code{"terms"}
#' @param na.action Not currently used (kept for a consistent signature with [predict.splm()])
#' @param ... Additional arguments (unused)
#'
#' @return The block-average prediction (the average of the point predictions
#'   over the rows of \code{newdata}), with a standard error and/or interval
#'   if requested
#'
#' @noRd
predict_block_splm <- function(object, newdata, se.fit, scale, df, interval, level, type, local, terms, na.action, ...) {
  # deal with local
  if (missing(local)) local <- NULL

  # check scale is numeric (if specified)
  if (!is.null(scale) && !is.numeric(scale)) {
    stop("scale must be numeric.", call. = FALSE)
  }

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
  if (object$dim_coords == 0) {
    # spcov_type "none"/"ie" fit without xcoord/ycoord: both coordinates are
    # synthetic placeholders (see get_point_ref_coords()), so newdata never
    # has them either -- fill both rather than requiring check_newdata_coords()
    # below to find columns that were never meant to exist in newdata
    obdata[[xcoord]] <- 0
    obdata[[ycoord]] <- 0
    newdata[[xcoord]] <- 0
    newdata[[ycoord]] <- 0
  } else if (object$dim_coords == 1) {
    obdata[[ycoord]] <- 0
    newdata[[ycoord]] <- 0
  }

  check_newdata_coords(newdata, xcoord, ycoord)

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

  newdata_model_list <- get_newdata_model_matrix(object, newdata)
  newdata <- newdata_model_list$newdata
  newdata_model <- newdata_model_list$newdata_model
  offset <- newdata_model_list$offset
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # finding rows w/out NA
  ob_predictors <- complete.cases(newdata_model)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }

  # a block prediction is the average of the point predictions over every row
  # of newdata, and because the model is linear in the fixed effects, that
  # average can be computed by first averaging the design matrix rows into a
  # single x0 (below) rather than predicting each point separately and
  # averaging afterward
  newdata_model_attr <- attributes(newdata_model)
  newdata_model <- matrix(colMeans(newdata_model), nrow = 1) # gives matrix, array class
  attr(newdata_model, "assign") <- newdata_model_attr$assign
  rownames(newdata_model) <- "1"
  colnames(newdata_model) <- newdata_model_attr$dimnames[[2]]
  x0 <- newdata_model
  betahat <- coef(object)
  cov_betahat <- vcov(object)
  y <- model.response(model.frame(object))
  offset <- model.offset(model.frame(object))
  # call terms if needed
  if (type == "terms") {
    return(predict_terms(object, newdata_model, se.fit, scale, df, interval, level, add_newdata_rows, terms, ...))
  }


  if (interval %in% c("none", "prediction")) {
    # local prediction list
    local <- get_local_list_prediction_block(local)

    # c0: each observed site's average covariance with the block (row-mean
    # of the newdata-rows-by-observed-sites covariance matrix), since the
    # block's covariance with an observed site is the average of that site's
    # covariance with every point in the block
    c0 <- colMeans(covmatrix(object, newdata = newdata, cov_type = "pred.obs"))
    Sig <- covmatrix(object)
    if (NROW(newdata) > 1e4) {
      # too many block points to form the dense NROW(newdata)^2 pred.pred
      # covariance matrix below, so compute the average pairwise covariance
      # row-by-row instead (see get_bk_var())
      s0 <- get_bk_var(object, newdata, local)
    } else {
      # s0: Var(block average) = average of all pairwise covariances among
      # block points (including each point's covariance with itself)
      s0 <- mean(covmatrix(object, newdata = newdata, cov_type = "pred.pred"))
    }
    Xmat <- model.matrix(object)

    if (local$method == "all") {
      cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(Sig)))
    } else {
      # restrict to the local$size observed sites closest (on average) to the
      # block, either by average distance or by average covariance with the
      # block, mirroring the point-prediction neighbor search in
      # get_pred_splm() but ranking against the whole block instead of a
      # single new location
      n <- length(c0)
      if (local$method == "distance") {
        dist_vector <- spdist_vectors(newdata, obdata, xcoord, ycoord, object$dim_coords)
        dist_vector <- colMeans(dist_vector)
        index <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, local$size))]
      } else if (local$method == "covariance") {
        # use abs() here for negative covariance types
        index <- order(abs(as.numeric(c0)))[seq(from = n, to = max(1, n - local$size + 1))]
      }
      obdata <- obdata[index, , drop = FALSE]
      c0 <- c0[index]
      Xmat <- Xmat[index, , drop = FALSE]
      y <- y[index]
      if (!is.null(offset)) {
        offset <- offset[index]
        y <- y - offset
      }
      cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(Sig[index, index, drop = FALSE])))
    }

    SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
    SqrtSigInv_y <- forwardsolve(cov_lowchol, y)
    residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
    SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)

    fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
    if (!is.null(offset)) {
      fit <- fit + offset
    }

    if (se.fit || interval == "prediction") {
      H <- x0 - crossprod(SqrtSigInv_c0, SqrtSigInv_X)
      vars <- as.numeric(s0 - crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% tcrossprod(cov_betahat, H))
      se <- sqrt(vars)
      if (!is.null(scale)) {
        se <- se * scale
        df <- df
      } else {
        df <- Inf
      }
      if (interval == "prediction") {
        tstar <- qt(1 - (1 - level) / 2, df = df)
        lwr <- fit - tstar * se
        upr <- fit + tstar * se
        fit <- cbind(fit, lwr, upr)
        row.names(fit) <- "1"
      }
      if (se.fit) {
        return(list(fit = fit, se.fit = se))
      } else {
        return(fit)
      }
    } else {
      return(fit)
    }
  } else if (interval == "confidence") {
    # finding fitted values of the mean parameters
    fit <- as.numeric(x0 %*% coef(object))
    # apply offset
    if (!is.null(offset)) {
      fit <- fit + offset
    }
    vars <- as.numeric(tcrossprod(x0 %*% cov_betahat, x0)) # different from
    # predict because x0 is a matrix here, not a vector
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
    row.names(fit) <- "1"
    if (se.fit) {
      return(list(fit = fit, se.fit = se))
    } else {
      return(fit)
    }
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}

#' Compute the block-average marginal variance for a large block (big data)
#'
#' @param object A fitted model object from [splm()]
#' @param newdata The data (points within the block) requiring prediction
#' @param local A fully-specified big-data \code{local} list
#'
#' @return The average pairwise covariance among the rows of \code{newdata}
#'   (i.e. \eqn{Var(\bar{Y})} for the block), computed row-by-row (optionally
#'   in parallel) rather than forming the full \code{NROW(newdata)}-by-\code{NROW(newdata)}
#'   covariance matrix, which would be too large for a big block
#'
#' @noRd
get_bk_var <- function(object, newdata, local) {
  index <- seq(1, NROW(newdata))
  object$obdata <- newdata
  if (local$parallel) {
    cl <- parallel::makeCluster(local$ncores)
    val <- parallel::parLapply(cl, index, get_each_bk_meancov, object, newdata)
    cl <- parallel::stopCluster(cl)
  } else {
    val <- lapply(index, get_each_bk_meancov, object, newdata)
  }
  mean(unlist(val))
}

#' Compute one row's average covariance with the rest of a prediction block
#'
#' @param index The row of \code{newdata} to compute the covariance row for
#' @param object A fitted model object from [splm()]
#' @param newdata The data (points within the block) requiring prediction
#'
#' @return The average of row \code{index}'s covariance with every row of
#'   \code{newdata} (including itself, with the independent error variance
#'   added back in, since \code{covmatrix()}'s \code{"obs.pred"} type omits it
#'   for a point predicted against itself)
#'
#' @noRd
get_each_bk_meancov <- function(index, object, newdata) {
  newdata <- newdata[index, , drop = FALSE]
  val <- as.vector(spmodel::covmatrix(object, newdata = newdata, cov_type = "obs.pred"))
  val[index] <- val[index] + object$coefficients$spcov[["ie"]] # this is to add ie variance to diagonal (which is omitted with newdata)
  mean(val)
}
