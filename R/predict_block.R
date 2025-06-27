predict_block_splm <- function(object, newdata, se.fit, scale, df, interval, level, type, local, terms, na.action, ...) {



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

  # finding rows w/out NA
  ob_predictors <- complete.cases(newdata_model)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }

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

    c0 <- colMeans(covmatrix(object, newdata = newdata, cov_type = "pred.obs"))
    Sig <- covmatrix(object)
    s0 <- mean(covmatrix(object, newdata = newdata, cov_type = "pred.pred"))
    Xmat <- model.matrix(object)

    if (local$method == "all") {
      cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(Sig)))
    } else {
      n <- length(c0)
      if (local$method == "distance") {
        # index <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, local$size))]
        index <- order(as.numeric(c0))[seq(from = n, to = max(1, n - local$size + 1))] # covariance method for now
      } else if (local$method == "covariance") {
        index <- order(as.numeric(c0))[seq(from = n, to = max(1, n - local$size + 1))]
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
