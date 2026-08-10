#' @param newdata_size The \code{size} value for each observation in \code{newdata}
#'   used when predicting for the binomial family, with a default value of 1.
#' @param var_correct A logical indicating whether to return the corrected prediction
#'   variances when predicting via models fit using \code{spglm()} or \code{spgautor()}. The default is
#'   \code{TRUE}.
#' @param delta A logical indicating whether to return delta method standard errors
#' on the response scale when \code{se.fit = TRUE} and \code{type = "response"}. The default is \code{FALSE}.
#' @param dispersion The dispersion of assumed when computing the prediction standard errors
#'   for \code{spglm()} or \code{spgautor()} model objects when \code{family}
#'   is \code{"nbinomial"}, \code{"beta"}, \code{"Gamma"}, or \code{"inverse.gaussian"}.
#'   If omitted, the model object dispersion parameter is used.
#' @rdname predict.spmodel
#' @method predict spglm
#' @order 9
#' @export
#' @examples
#' \donttest{
#' spgmod <- spglm(presence ~ elev * strat,
#'   family = "binomial",
#'   data = moose,
#'   spcov_type = "exponential"
#' )
#' predict(spgmod, moose_preds)
#' predict(spgmod, moose_preds, interval = "prediction")
#' augment(spgmod, newdata = moose_preds, interval = "prediction")
#' }
predict.spglm <- function(object, newdata, type = c("link", "response", "terms", "weight"), se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                          level = 0.95, dispersion = NULL, terms = NULL, local, var_correct = TRUE, delta = FALSE, newdata_size, na.action = na.fail, ...) {
  # match type argument so the two display
  type <- match.arg(type)
  if (type == "weight") {
    se.fit <- FALSE
    interval <- "none"
  }

  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with newdata_size
  if (missing(newdata_size)) newdata_size <- NULL

  # deal with local
  if (missing(local)) local <- NULL

  if (!is.logical(delta)) {
    stop("delta must be TRUE or FALSE", call. = FALSE) # consider making delta relevant default to match glm
  }

  # build the prediction setup object and pull its elements into named local
  # variables (explicit assignment, not list2env(), so static analysis --
  # R CMD check's codetools-based check and RStudio's diagnostics -- can see
  # where object/obdata/xcoord/ycoord/newdata/etc. below come from)
  prediction_object <- get_prediction_object_spglm(object, newdata, dispersion, newdata_size, local)
  object <- prediction_object$object
  local <- prediction_object$local
  newdata_size <- prediction_object$newdata_size
  obdata <- prediction_object$obdata
  xcoord <- prediction_object$xcoord
  ycoord <- prediction_object$ycoord
  newdata <- prediction_object$newdata
  add_newdata_rows <- prediction_object$add_newdata_rows
  spcov_params_val <- prediction_object$spcov_params_val
  dispersion_params_val <- prediction_object$dispersion_params_val
  randcov_params_val <- prediction_object$randcov_params_val
  newdata_model <- prediction_object$newdata_model
  offset <- prediction_object$offset

  # call terms if needed
  if (type == "terms") {
    # glm supports standard errors for terms objects but not intervals (no interval argument)
    # scale df not used for glms
    return(predict_terms(object, newdata_model, se.fit, scale = NULL, df = Inf, interval, level, add_newdata_rows, terms, ...))
  }

  # storing newdata as a list
  npred <- NROW(newdata)
  newdata_rows_list <- split(newdata, seq_len(npred))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(npred))

  # storing newdata as a list (row_index lets get_pred_spglm() slice its row
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

    # extend the prediction object with everything get_pred_spglm() needs that
    # is constant across prediction rows, so it can rely on a single object
    # argument instead of two dozen individually named ones. assignment must
    # go through `[names(.)] <-` rather than `$<-`/`[[<-`,
    # since those delete a key entirely when its value is NULL (e.g. random,
    # partition_factor, cov_lowchol) instead of storing the NULL
    pred_row_context <- list(
      se.fit = se.fit, interval = interval, formula = object$terms,
      random = object$random, randcov_terms = randcov_terms,
      partition_factor = object$partition_factor, reform_bar2 = reform_bar2,
      partition_index_obdata = partition_index_obdata, cov_lowchol = cov_lowchol,
      Xmat = model.matrix(object), y = object$y, dim_coords = object$dim_coords,
      betahat = coefficients(object), cov_betahat = vcov(object, var_correct = FALSE),
      contrasts = object$contrasts, local = local_list, family = object$family,
      w = fitted(object, type = "link"), model_offset = model.offset(model.frame(object)),
      size = object$size, dispersion = dispersion_params_val,
      predvar_adjust_ind = predvar_adjust_ind, xlevels = object$xlevels,
      diagtol = object$diagtol, type = type,
      dist_matrix_full = dist_matrix_full, partition_vector_full = partition_vector_full,
      cov_vector_full = cov_vector_full
    )
    prediction_object[names(pred_row_context)] <- pred_row_context

    pred_spglm <- run_pred_dispatch(get_pred_spglm, newdata_list, prediction_object, local_list)

    if (type == "weight") {
      fit <- do.call("rbind", lapply(pred_spglm, function(x) x$fit))
      if (add_newdata_rows) {
        colnames(fit) <- object$observed_index
        rownames(fit) <- object$missing_index
      }
      return(fit)
    }

    if (interval == "none") {
      fit <- vapply(pred_spglm, function(x) x$fit, numeric(1))
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
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
        if (type == "response" && se.fit && delta) {
          se <- get_delta_se(fit, se, object$family, newdata_size)
        }
        if (type == "response") {
          fit <- invlink(fit, object$family, newdata_size)
        }
      } else {
        se <- NULL
        if (type == "response") {
          fit <- invlink(fit, object$family, newdata_size)
        }
      }
      return(finalize_interval_none(fit, se, add_newdata_rows, object$missing_index))
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
      tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      if (type == "response" && se.fit && delta) {
        se <- get_delta_se(fit, se, object$family, newdata_size)
      }
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
        lwr <- invlink(lwr, object$family, newdata_size)
        upr <- invlink(upr, object$family, newdata_size)
      }
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
    tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    if (type == "response") {
      fit <- invlink(fit, object$family, newdata_size)
      lwr <- invlink(lwr, object$family, newdata_size)
      upr <- invlink(upr, object$family, newdata_size)
    }
    return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows, object$missing_index))
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}

#' Predict a single new observation for an \code{spglm()} model
#'
#' @param newdata_list A list with elements \code{row} (a single-row data frame
#'   for the new observation), \code{x0} (its design matrix row), and
#'   \code{row_index} (its row number within \code{newdata}, used to slice
#'   \code{dist_matrix_full}/\code{partition_vector_full}/\code{cov_vector_full}
#'   when those are supplied)
#' @param prediction_object A list of values constant across prediction rows
#'   (built by \code{predict.spglm()}), assigned into named local variables at
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
#'   \code{y} (the observed response vector, unused directly since \code{w}
#'   carries the fitted latent values), \code{betahat}, \code{cov_betahat},
#'   \code{dim_coords}, \code{contrasts}, \code{local} (a fully-specified
#'   big-data \code{local} list), \code{family}, \code{w} (the observed-data
#'   latent link-scale predictor vector, with the offset still added),
#'   \code{model_offset} (the observed-data offset, or \code{NULL}), \code{size}
#'   (binomial trial sizes, used only when \code{family} is \code{"binomial"}),
#'   \code{dispersion}, \code{predvar_adjust_ind} (whether to add the
#'   \code{get_wts_varw()} adjustment for estimation of the latent random
#'   effects to the prediction variance), \code{xlevels}, \code{diagtol},
#'   \code{type} (\code{"link"}, \code{"response"}, or \code{"weight"}), and,
#'   for medium-sized data (see \code{predict.spglm()}), \code{dist_matrix_full},
#'   \code{partition_vector_full}, and \code{cov_vector_full} -- the
#'   observed-by-prediction distance, partition, and covariance matrices,
#'   precomputed once (vectorized across every prediction row) so this row's
#'   values can be sliced out of them instead of reconstructed from scratch;
#'   \code{NULL} for big data, where reconstructing them one row at a time
#'   instead of all at once avoids ever holding an observed-by-prediction
#'   matrix in memory
#'
#' @return A list with element \code{fit} (and, if \code{se.fit} or
#'   \code{interval == "prediction"}, element \code{var}) for the new
#'   observation, on the link scale
#'
#' @noRd
get_pred_spglm <- function(newdata_list, prediction_object) {
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
  betahat <- prediction_object$betahat
  cov_betahat <- prediction_object$cov_betahat
  dim_coords <- prediction_object$dim_coords
  contrasts <- prediction_object$contrasts
  local <- prediction_object$local
  family <- prediction_object$family
  w <- prediction_object$w
  model_offset <- prediction_object$model_offset
  size <- prediction_object$size
  dispersion <- prediction_object$dispersion
  predvar_adjust_ind <- prediction_object$predvar_adjust_ind
  xlevels <- prediction_object$xlevels
  diagtol <- prediction_object$diagtol
  type <- prediction_object$type
  dist_matrix_full <- prediction_object$dist_matrix_full
  partition_vector_full <- prediction_object$partition_vector_full
  cov_vector_full <- prediction_object$cov_vector_full

  # adjust w: strip the offset back out of the fitted latent link-scale
  # values so the covariance-based prediction below operates on the same
  # offset-free scale as the whitened X/c0 (the offset for this new
  # observation is added back in by the caller, e.g. predict.spglm())
  if (!is.null(model_offset)) {
    w <- w - model_offset
  }


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
    w <- w[nn_index]
    y <- y[nn_index]
    if (!is.null(size)) {
      size <- size[nn_index]
    }
  }


  if (local$method == "covariance") {
    n <- length(cov_vector_val)
    # use abs here for the most covariance
    cov_index <- order(abs(as.numeric(cov_vector_val)))[seq(from = n, to = max(1, n - local$size + 1))]
    obdata <- obdata[cov_index, , drop = FALSE]
    cov_vector_val <- cov_vector_val[cov_index]
    w <- w[cov_index]
    y <- y[cov_index]
    if (!is.null(size)) {
      size <- size[cov_index]
    }
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
  }

  # "whiten" X, w, and c0 by left-multiplying by Sigma^{-1/2} -- implemented
  # as forward substitution against the lower Cholesky factor rather than
  # explicitly forming Sigma^{-1}
  c0 <- as.numeric(cov_vector_val)
  SqrtSigInv_X <- base::forwardsolve(cov_lowchol, Xmat)
  SqrtSigInv_w <- base::forwardsolve(cov_lowchol, w)
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
    # universal kriging BLUP on the link scale: the trend x0 %*% betahat plus
    # a covariance-weighted combination of the observed (whitened) latent
    # residuals
    residuals_pearson <- SqrtSigInv_w - SqrtSigInv_X %*% betahat
    fit <- as.numeric(x0 %*% betahat + base::crossprod(SqrtSigInv_c0, residuals_pearson))
  }

  if (se.fit || interval == "prediction") {
    # kriging prediction variance (see the analogous comment in predict.R's
    # get_pred_splm()): marginal variance, minus variance explained by the
    # observed data, plus the betahat-estimation-uncertainty inflation term
    H <- x0 - base::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
    # a random slope's contribution to Var(Y0) is sigma^2 * x0^2, not sigma^2
    # (as it would be for a random intercept), so it must be computed for this
    # specific newdata row rather than summed directly from randcov_params_val
    total_var <- spcov_params_val[["de"]] + spcov_params_val[["ie"]] +
      randcov_newvar(randcov_params_val, newdata_list$row, randcov_terms)
    var <- as.numeric(total_var - base::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% base::tcrossprod(cov_betahat, H))
    if (predvar_adjust_ind) {
      # for GLMs, w itself is a latent quantity estimated (not observed
      # directly like y in the Gaussian case), so its own estimation
      # uncertainty contributes an extra variance term on top of the usual
      # kriging variance above
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
predict.spgautor <- function(object, newdata, type = c("link", "response", "terms", "weight"), se.fit = FALSE,
                             interval = c("none", "confidence", "prediction"),
                             level = 0.95, dispersion = NULL, terms = NULL, local, var_correct = TRUE, delta = FALSE, newdata_size, na.action = na.fail, ...) {
  # match type argument so the two display
  type <- match.arg(type)
  if (type == "weight") {
    se.fit <- FALSE
    interval <- "none"
  }

  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with newdata_size
  if (missing(newdata_size)) newdata_size <- NULL

  # deal with local
  if (missing(local)) local <- NULL

  if (!is.logical(delta)) {
    stop("delta must be TRUE or FALSE", call. = FALSE) # consider making delta relevant default to match glm
  }

  # build the prediction setup object and pull its elements into named local
  # variables (explicit assignment, not list2env(), so static analysis --
  # R CMD check's codetools-based check and RStudio's diagnostics -- can see
  # where object/newdata/spcov_params_val/newdata_model/etc. below come from)
  prediction_object <- get_prediction_object_spgautor(object, newdata, dispersion, newdata_size, local)
  object <- prediction_object$object
  local <- prediction_object$local
  newdata_size <- prediction_object$newdata_size
  newdata <- prediction_object$newdata
  spcov_params_val <- prediction_object$spcov_params_val
  dispersion_params_val <- prediction_object$dispersion_params_val
  randcov_params_val <- prediction_object$randcov_params_val
  newdata_model <- prediction_object$newdata_model
  offset <- prediction_object$offset

  # call terms if needed
  if (type == "terms") {
    # scale df not used for glms
    return(predict_terms(object, newdata_model, se.fit, scale = NULL, df = Inf, interval, level, add_newdata_rows = TRUE, terms, ...))
  }


  # storing newdata as a list
  newdata_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  if (interval %in% c("none", "prediction")) {
    # like predict.spautor(), spgautor() models don't use a big-data local
    # approximation here: the fixed neighbor structure (object$W) already
    # gives one joint covariance matrix over all sites, built once and then
    # sliced by row/column index below
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

    # lower triangular cholesky (dense conversion avoids Matrix S4 chol()
    # returning a factorization object that base::t() in get_wts_varw() cannot transpose)
    cov_matrix_lowchol <- t(base::chol(as.matrix(cov_matrix_val[object$observed_index, object$observed_index, drop = FALSE])))

    # find X observed
    X <- model.matrix(object)
    SqrtSigInv_X <- forwardsolve(cov_matrix_lowchol, X)

    # find w observed
    w <- fitted(object, type = "link")
    model_offset <- model.offset(model.frame(object))
    if (!is.null(model_offset)) {
      w <- w - model_offset
    }
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

    # storing each new observation's row-specific quantities as a list
    cluster_list <- mapply(
      x0 = newdata_model_list, c0 = cov_vector_val_list, s0 = total_var_list,
      FUN = function(x0, c0, s0) list(x0 = x0, c0 = c0, s0 = s0), SIMPLIFY = FALSE
    )

    # extend the prediction object with everything get_pred_spgautor() needs
    # that is constant across prediction rows, so it can rely on a single
    # object argument instead of many individually named ones
    pred_row_context <- list(
      cov_matrix_lowchol = cov_matrix_lowchol, betahat = betahat,
      residuals_pearson_w = residuals_pearson_w, cov_betahat = cov_betahat,
      SqrtSigInv_X = SqrtSigInv_X, se.fit = se.fit, interval = interval,
      type = type, Xmat = X
    )
    prediction_object[names(pred_row_context)] <- pred_row_context

    pred_spautor <- run_pred_dispatch(get_pred_spgautor, cluster_list, prediction_object, local_list)

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
        if (type == "response" && se.fit && delta) {
          se <- get_delta_se(fit, se, object$family, newdata_size)
        }
        if (type == "response") {
          fit <- invlink(fit, object$family, newdata_size)
        }
      } else {
        se <- NULL
        if (type == "response") {
          fit <- invlink(fit, object$family, newdata_size)
        }
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
      tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      if (type == "response" && se.fit && delta) {
        se <- get_delta_se(fit, se, object$family, newdata_size)
      }
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
        lwr <- invlink(lwr, object$family, newdata_size)
        upr <- invlink(upr, object$family, newdata_size)
      }
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
    tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    if (type == "response") {
      fit <- invlink(fit, object$family, newdata_size)
      lwr <- invlink(lwr, object$family, newdata_size)
      upr <- invlink(upr, object$family, newdata_size)
    }
    return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows = TRUE, object$missing_index))
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}

#' Predict a single new observation for an \code{spgautor()} model
#'
#' @param cluster_list A list with elements \code{x0} (the design matrix row
#'   for the new observation), \code{c0} (the covariance vector between the
#'   new observation and the observed data), and \code{s0} (the marginal
#'   variance of the new observation)
#' @param prediction_object A list of values constant across prediction rows
#'   (built by \code{predict.spgautor()}), assigned into named local
#'   variables at the top of this function: \code{cov_matrix_lowchol} (the lower Cholesky factor of
#'   the observed-data covariance matrix), \code{betahat}, \code{residuals_pearson_w}
#'   (Pearson residuals on the link scale for the observed data), \code{cov_betahat},
#'   \code{SqrtSigInv_X} (\code{cov_matrix_lowchol^{-1} \%*\% Xmat}), \code{se.fit},
#'   \code{interval} (\code{"none"}, \code{"confidence"}, or \code{"prediction"}),
#'   \code{type} (\code{"response"} or \code{"weight"}), and \code{Xmat} (the
#'   observed-data design matrix)
#'
#' @return A list with element \code{fit} (and, if \code{se.fit} or
#'   \code{interval == "prediction"}, element \code{var}) for the new
#'   observation, on the link scale
#'
#' @noRd
get_pred_spgautor <- function(cluster_list, prediction_object) {
  # explicit assignment, not list2env(), so static analysis (R CMD check's
  # codetools-based check and RStudio's diagnostics) can see where each name
  # below comes from
  cov_matrix_lowchol <- prediction_object$cov_matrix_lowchol
  betahat <- prediction_object$betahat
  residuals_pearson_w <- prediction_object$residuals_pearson_w
  cov_betahat <- prediction_object$cov_betahat
  SqrtSigInv_X <- prediction_object$SqrtSigInv_X
  se.fit <- prediction_object$se.fit
  interval <- prediction_object$interval
  type <- prediction_object$type
  Xmat <- prediction_object$Xmat

  # shared with get_pred_spautor() -- see get_areal_pred() in
  # R/predict_helpers.R. Note: no get_wts_varw() adjustment here, unlike
  # get_pred_spglm(), since spgautor's var_correct adjustment is instead
  # applied by the caller, predict.spgautor(), after this function returns
  get_areal_pred(cluster_list, cov_matrix_lowchol, betahat, residuals_pearson_w, cov_betahat, SqrtSigInv_X, se.fit, interval, type, Xmat)
}

#' @name predict.spmodel
#' @method predict spglm_list
#' @order 11
#' @export
predict.spglm_list <- function(object, newdata, type = c("link", "response", "terms"), se.fit = FALSE,
                               interval = c("none", "confidence", "prediction"), level = 0.95,
                               dispersion = NULL, terms = NULL, local, var_correct = TRUE, newdata_size,
                               na.action = na.fail, ...) {
  type <- match.arg(type)
  # match interval argument so the three display
  interval <- match.arg(interval)

  # deal with local
  if (missing(local)) local <- NULL

  # deal with newdata_size
  if (missing(newdata_size)) newdata_size <- NULL

  if (missing(newdata)) {
    preds <- lapply(object, function(x) {
      predict(
        x,
        type = type,
        se.fit = se.fit,
        interval = interval,
        level = level,
        dispersion = dispersion,
        terms = terms,
        local = local,
        var_correct = var_correct,
        newdata_size = newdata_size, # don't need na.action as it is fixed in predict
        ...
      )
    })
  } else {
    preds <- lapply(object, function(x) {
      predict(
        x,
        newdata = newdata,
        type = type,
        se.fit = se.fit,
        interval = interval,
        level = level,
        dispersion = dispersion,
        terms = terms,
        local = local,
        var_correct = var_correct,
        newdata_size = newdata_size, # don't need na.action as it is fixed in predict
        ...
      )
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
