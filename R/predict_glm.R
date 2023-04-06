#' @param type The scale (\code{response} or \code{link}) of predictions obtained
#'   using \code{spglm()} or \code{spgautor} objects.
#' @param newdata_size The \code{size} value for each observation in \code{newdata}
#'   used when predicting for the binomial family.
#' @param var_correct A logical indicating whether to return the corrected prediction
#'   variances when predicting via models fit using \code{spglm()} or \code{spgautor()}. The default is
#'   \code{TRUE}.
#' @rdname predict.spmodel
#' @method predict spglm
#' @order 9
#' @export
predict.spglm <- function(object, newdata, type = c("link", "response"), se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                         newdata_size, level = 0.95, local, var_correct = TRUE, ...) {


  # match type argument so the two display
  type <- match.arg(type)

  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with newdata_size
  if (missing(newdata_size)) newdata_size <- NULL

  # deal with local
  if (missing(local)) {
    local <- NULL
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

  # set newdata_size if needed
  if (is.null(newdata_size) && object$family == "binomial") {
    newdata_size <- rep(1, NROW(newdata))
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

  # save dispersion param vector
  dispersion_params_val <- as.vector(coef(object, type = "dispersion")) # remove class

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

  # storing newdata as a list
  newdata_rows_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_list <- mapply(x = newdata_rows_list, y = newdata_model_list, FUN = function(x, y) list(row = x, x0 = y), SIMPLIFY = FALSE)

  if (interval %in% c("none", "prediction")) {



    # local prediction list
    local_list <- get_local_list_prediction(local)

    # random stuff
    if (!is.null(object$random)) {
      randcov_names <- get_randcov_names(object$random)
      randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
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
      partition_index_obdata <- list(reform_bar2_vals = p_reform_bar2_vals, reform_bar2_xlev = p_reform_bar2_xlev)
      # partition_index_obdata <- as.vector(model.matrix(reform_bar2, obdata))
    } else {
      reform_bar2 <- NULL
      partition_index_obdata <- NULL
    }

    # matrix cholesky
    if (local_list$method == "all") {
      cov_matrix_val <- covmatrix(object)
      cov_lowchol <- t(chol(cov_matrix_val))
      predvar_adjust_ind <- FALSE
      predvar_adjust_all <- TRUE
    } else {
      cov_lowchol <- NULL
      predvar_adjust_ind <- TRUE
      predvar_adjust_all <- FALSE
    }

    # change predvar adjust based on var correct
    if (!var_correct) {
      predvar_adjust_ind <- FALSE
      predvar_adjust_all <- FALSE
    }

    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      pred_spglm <- parallel::parLapply(cl, newdata_list, get_pred_spglm,
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
                                       y = object$y, dim_coords = object$dim_coords,
                                       betahat = coefficients(object), cov_betahat = vcov(object, var_correct = FALSE),
                                       contrasts = object$contrasts,
                                       local = local_list, family = object$family, w = fitted(object, type = "link"), size = object$size,
                                       dispersion = dispersion_params_val, predvar_adjust_ind = predvar_adjust_ind, diagtol = object$diagtol
      )
      cl <- parallel::stopCluster(cl)
    } else {
      pred_spglm <- lapply(newdata_list, get_pred_spglm,
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
                          y = object$y, dim_coords = object$dim_coords,
                          betahat = coefficients(object), cov_betahat = vcov(object, var_correct = FALSE),
                          contrasts = object$contrasts,
                          local = local_list, family = object$family,
                          w = fitted(object, type = "link"), size = object$size,
                          dispersion = dispersion_params_val, predvar_adjust_ind = predvar_adjust_ind,
                          diagtol = object$diagtol
      )
    }

    if (interval == "none") {
      fit <- vapply(pred_spglm, function(x) x$fit, numeric(1))
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
      }
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
      }
      if (se.fit) {
        vars <- vapply(pred_spglm, function(x) x$var, numeric(1))
        if (predvar_adjust_all) {
          # predvar_adjust is for the local function so FALSE there is TRUE
          # here
          vars_adj <- get_wts_varw(
            family = object$family,
            Xmat = model.matrix(object),
            y = object$y,
            w = fitted(object, type = "link"),
            size = object$size,
            dispersion = dispersion_params_val,
            cov_lowchol = cov_lowchol,
            x0 = newdata_model,
            c0 = covmatrix(object, newdata)
          )
          vars <- vars_adj + vars
        }
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
      fit <- vapply(pred_spglm, function(x) x$fit, numeric(1))
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
      }
      vars <- vapply(pred_spglm, function(x) x$var, numeric(1))
      if (predvar_adjust_all) {
        vars_adj <- get_wts_varw(
          family = object$family,
          Xmat = model.matrix(object),
          y = object$y,
          w = fitted(object, type = "link"),
          size = object$size,
          dispersion = dispersion_params_val,
          cov_lowchol = cov_lowchol,
          x0 = newdata_model,
          c0 = covmatrix(object, newdata)
        )
        vars <- vars_adj + vars
      }
      se <- sqrt(vars)
      # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
      tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
        lwr <- invlink(lwr, object$family, newdata_size)
        upr <- invlink(upr, object$family, newdata_size)
      }
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
    # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
    tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    if (type == "response") {
      fit <- invlink(fit, object$family, newdata_size)
      lwr <- invlink(lwr, object$family, newdata_size)
      upr <- invlink(upr, object$family, newdata_size)
    }
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

get_pred_spglm <- function(newdata_list, se.fit, interval, formula, obdata, xcoord, ycoord,
                          spcov_params_val, random, randcov_params_val, reform_bar2_list,
                          Z_index_obdata_list, reform_bar1_list, Z_val_obdata_list, partition_factor,
                          reform_bar2, partition_index_obdata, cov_lowchol,
                          Xmat, y, betahat, cov_betahat, dim_coords, contrasts, local,
                          family, w, size, dispersion, predvar_adjust_ind, diagtol) {


  # storing partition vector
  partition_vector <- partition_vector(partition_factor,
                                       data = obdata,
                                       newdata = newdata_list$row, reform_bar2 = reform_bar2,
                                       partition_index_data = partition_index_obdata
  )

  # subsetting partition vector
  if (!is.null(partition_vector) && local$method %in% c("distance", "covariance") &&
      !labels(terms(partition_factor)) %in% labels(terms(random))) {
    obdata <- obdata[as.vector(partition_vector) == 1, , drop = FALSE]
    partition_vector <- Matrix(1, nrow = 1, ncol = NROW(obdata))
  }

  dist_vector <- spdist_vectors(newdata_list$row, obdata, xcoord, ycoord, dim_coords)

  # subsetting data if method distance
  if (local$method == "distance") {
    n <- length(dist_vector)
    nn_index <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, local$size))]
    obdata <- obdata[nn_index, , drop = FALSE]
    dist_vector <- dist_vector[, nn_index]
    w <- w[nn_index]
    if (!is.null(size)) {
      size <- size[nn_index]
    }
  }

  # making random vector if necessary
  if (!is.null(randcov_params_val)) {
    randcov_vector_val <- randcov_vector(randcov_params_val, obdata, newdata_list$row, reform_bar2_list, Z_index_obdata_list)
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
    w <- w[cov_index]
    if (!is.null(size)) {
      size <- size[cov_index]
    }
  }

  if (local$method %in% c("distance", "covariance")) {
    if (!is.null(random)) {
      randcov_names <- get_randcov_names(random)
      randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
    }
    partition_matrix_val <- partition_matrix(partition_factor, obdata)
    cov_matrix_val <- cov_matrix(
      spcov_params_val, spdist(obdata, xcoord, ycoord), randcov_params_val,
      randcov_Zs, partition_matrix_val, diagtol = diagtol
    )
    cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(cov_matrix_val)))
    model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass)
    Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
    y <- model.response(model_frame)
    if (NCOL(y) == 2) {
      y_modr <- y
      y <- y_modr[, 1, drop = FALSE]
      size <- rowSums(y_modr)
    } else {
      if (family == "binomial") {
        size <- rep(1, n)
      } else {
        size <- NULL
      }
    }
  }



  c0 <- as.numeric(cov_vector_val)
  SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
  SqrtSigInv_w <- forwardsolve(cov_lowchol, w)
  residuals_pearson <- SqrtSigInv_w - SqrtSigInv_X %*% betahat
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
  x0 <- newdata_list$x0

  fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
  H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  if (se.fit || interval == "prediction") {
    total_var <- sum(spcov_params_val[["de"]], spcov_params_val[["ie"]], randcov_params_val)
    var <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
    if (predvar_adjust_ind) {
      var_adj <- get_wts_varw(family, Xmat, y, w, size, dispersion, cov_lowchol, x0, c0)
      var <- var_adj + var
    }
    pred_list <- list(fit = fit, var = var)
  } else {
    pred_list <- list(fit = fit)
  }
  pred_list
}







#' @rdname predict.spmodel
#' @method predict spgautor
#' @order 10
#' @export
predict.spgautor <- function(object, newdata, type = c("link", "response"), se.fit = FALSE,
                             interval = c("none", "confidence", "prediction"),
                            newdata_size, level = 0.95, local, var_correct = TRUE, ...) {

  # match type argument so the two display
  type <- match.arg(type)

  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with newdata_size
  if (missing(newdata_size)) newdata_size <- NULL

  # deal with local
  if (missing(local)) {
    local <- NULL
  }

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }

  # deal with local
  if (is.null(local)) {
    local <- FALSE
  }

  # write newdata if predicting missing data
  newdata <- object$data[object$missing_index, , drop = FALSE]

  # set newdata_size if needed
  if (is.null(newdata_size) && object$family == "binomial") {
    newdata_size <- rep(1, NROW(newdata))
  }

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save dispersion param vector
  dispersion_params_val <- as.vector(coef(object, type = "dispersion")) # remove class

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

    # # lower triangular cholesky
    cov_matrix_lowchol <- t(chol(cov_matrix_val[object$observed_index, object$observed_index, drop = FALSE]))
    # cov_matrix_lowchol <- t(chol(cov_matrix_val_obs))

    # find X observed
    X <- model.matrix(object)
    SqrtSigInv_X <- forwardsolve(cov_matrix_lowchol, X)

    # find w observed
    w <- fitted(object, type = "link")
    SqrtSigInv_w <- forwardsolve(cov_matrix_lowchol, w)

    # beta hat
    betahat <- coef(object)

    # residuals pearson
    residuals_pearson_w <- SqrtSigInv_w - SqrtSigInv_X %*% betahat

    # cov beta hat
    cov_betahat <- vcov(object, var_correct = FALSE)

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
      pred_spautor <- parallel::parLapply(cl, cluster_list, get_pred_spgautor_parallel,
                                          cov_matrix_lowchol, betahat,
                                          residuals_pearson_w,
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
          get_pred_spgautor(
            x0 = x0, c0 = c0, s0 = s0,
            cov_matrix_lowchol, betahat,
            residuals_pearson_w,
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
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
      }
      if (se.fit) {
        vars <- vapply(pred_spautor, function(x) x$var, numeric(1))
        if (var_correct) {
          vars_adj <- get_wts_varw(
            family = object$family,
            Xmat = model.matrix(object),
            y = object$y,
            w = fitted(object, type = "link"),
            size = object$size,
            dispersion = dispersion_params_val,
            cov_lowchol = cov_matrix_lowchol,
            x0 = newdata_model,
            c0 = cov_vector_val
          )
          vars <- vars_adj + vars
        }
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
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
      }
      vars <- vapply(pred_spautor, function(x) x$var, numeric(1))
      if (var_correct) {
        vars_adj <- get_wts_varw(
          family = object$family,
          Xmat = model.matrix(object),
          y = object$y,
          w = fitted(object, type = "link"),
          size = object$size,
          dispersion = dispersion_params_val,
          cov_lowchol = cov_matrix_lowchol,
          x0 = newdata_model,
          c0 = cov_vector_val
        )
        vars <- vars_adj + vars
      }
      se <- sqrt(vars)
      # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
      tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
        lwr <- invlink(lwr, object$family, newdata_size)
        upr <- invlink(upr, object$family, newdata_size)
      }
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
    # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
    tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    if (type == "response") {
      fit <- invlink(fit, object$family, newdata_size)
      lwr <- invlink(lwr, object$family, newdata_size)
      upr <- invlink(upr, object$family, newdata_size)
    }
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

get_pred_spgautor <- function(x0, c0, s0, cov_matrix_lowchol, betahat, residuals_pearson_w, cov_betahat, SqrtSigInv_X, se.fit, interval) {
  SqrtSigInv_c0 <- forwardsolve(cov_matrix_lowchol, c0)
  fit <- as.numeric(x0 %*% betahat + crossprod(SqrtSigInv_c0, residuals_pearson_w))
  if (se.fit || interval == "prediction") {
    H <- x0 - crossprod(SqrtSigInv_c0, SqrtSigInv_X)
    var <- as.numeric(s0 - crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% tcrossprod(cov_betahat, H))
    pred_list <- list(fit = fit, var = var)
  } else {
    pred_list <- list(fit = fit)
  }
  pred_list
}

get_pred_spgautor_parallel <- function(cluster_list, cov_matrix_lowchol, betahat, residuals_pearson_w, cov_betahat, SqrtSigInv_X, se.fit, interval) {
  x0 <- cluster_list$x0
  c0 <- cluster_list$c0
  s0 <- cluster_list$s0
  get_pred_spgautor(x0, c0, s0, cov_matrix_lowchol, betahat, residuals_pearson_w, cov_betahat, SqrtSigInv_X, se.fit, interval)
}

#' @name predict.spmodel
#' @method predict spglm_list
#' @order 11
#' @export
predict.spglm_list <- function(object, newdata, type = c("link", "response"), se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                              newdata_size, level = 0.95, local, var_correct = TRUE, ...) {

  type <- match.arg(type)
  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with local
  if (missing(local)) {
    local <- NULL
  }

  # deal with newdata_size
  if (missing(newdata_size)) {
    newdata_size <- NULL
  }

  if (missing(newdata)) {
    preds <- lapply(object, function(x) {
      predict(x, type = type, se.fit = se.fit, interval = interval, newdata_size = newdata_size, level = level, local = local, var_correct = var_correct, ...)
    })
  } else {
    preds <- lapply(object, function(x) {
      predict(x, newdata = newdata, type = type, se.fit = se.fit, interval = interval, newdata_size = newdata_size, level = level, local = local, var_correct = var_correct, ...)
    })
  }
  names(preds) <- names(object)
  preds
}

#' @name predict.spmodel
#' @method predict spgautor_list
#' @order 12
#' @export
predict.spgautor_list <- predict.spglm_list
