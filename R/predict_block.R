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

  # deal with local, which will be additionally processed further down
  local_unset <- is.null(local)
  if (local_unset) local <- FALSE

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
  offset_newdata <- newdata_model_list$offset
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
  model_frame <- model.frame(object)
  y <- model.response(model_frame)
  # Two different offsets are in play and should not not be confused. The
  # observed-data offset comes out of y so the kriging below runs on the
  # offset-free scale. The block's own offset is the average of the newdata offsets,
  # because a block prediction is the average of the point predictions and each
  # point prediction carries its own offset; it goes back on at the end,
  # alongside the same colMeans() averaging already applied to newdata_model.
  offset_obdata <- model.offset(model_frame)
  if (!is.null(offset_obdata)) {
    y <- y - as.vector(offset_obdata)
  }
  offset_block <- if (is.null(offset_newdata)) NULL else mean(as.vector(offset_newdata))
  # call terms if needed
  if (type == "terms") {
    return(predict_terms(object, newdata_model, se.fit, scale, df, interval, level, add_newdata_rows, terms, ...))
  }


  if (interval %in% c("none", "prediction")) {
    # finish the unset-local defaulting
    if (local_unset) {
      approx_obs <- object$n > 10000
      approx_pred <- NROW(newdata) > 10000
      if (approx_obs || approx_pred) {
        local <- list(
          method = if (approx_obs) "covariance" else "all",
          size = 4000,
          method_new = "basis",
          size_new = if (approx_pred) 4000L else Inf,
          ordering = "grts"
        )
        sides <- c(
          if (approx_obs) "the fitted model sample size",
          if (approx_pred) "the number of prediction locations in newdata"
        )
        message(sprintf(
          paste0(
            "Because %s %s 10,000, we are using a computationally efficient ",
            "approximation for the block prediction. To compute the exact solution instead, rerun ",
            "predict() with local = FALSE. Be aware that local = FALSE may result in ",
            "exceedingly long computational times."
          ),
          paste(sides, collapse = " and "),
          if (length(sides) > 1L) "each exceed" else "exceeds"
        ))
      }
    }

    # local prediction list
    local <- get_local_list_prediction_block(local)

    Sig <- covmatrix(object)

    # c0: each observed site's average covariance with the block (column-mean
    # of the newdata-rows-by-observed-sites covariance matrix), since the
    # block's covariance with an observed site is the average of that site's
    # covariance with every point in the block.
    # s0: Var(block average) = the average pairwise covariance among the block
    # points (each point's covariance with itself included).
    G <- NROW(newdata)
    if (local$size_new >= G) {
      # local = FALSE, or size_new >= G: every block point is an s0 node, so
      # c0, x0, and s0 are all exact
      bq <- get_block_quantities(object, newdata,
        nodes = seq_len(G),
        parallel = local$parallel, ncores = local$ncores
      )
      c0 <- bq$c0
      s0 <- bq$s0
    } else {
      # size_new well-spread nodes chosen from newdata via the same
      # get_decorrelate_order() helper decorrelate()/sprnorm() use;
      # coordinates here are already anisotropy-transformed. The default is
      # "grts".
      if (object$dim_coords == 0) {
        # no meaningful coordinates ("none"/"ie" covariance)
        nodes <- seq_len(local$size_new)
      } else {
        block_ordering <- local$ordering
        if (is.null(block_ordering)) {
          block_ordering <- "grts"
        }
        block_ord <- get_decorrelate_order(block_ordering, newdata[[xcoord]], newdata[[ycoord]])$order
        nodes <- block_ord[seq_len(local$size_new)]
      }

      if (identical(local$method_new, "subset")) {
        # x0 (computed over the full newdata) stays exact; c0 and s0 are
        # taken on the size_new nodes only; s0's diagonal lands at weight
        # 1 / size_new here; re-weight it to the 1 / G it carries in the exact
        # block variance.
        bq <- get_block_quantities(object, newdata[nodes, , drop = FALSE],
          nodes = seq_len(local$size_new),
          parallel = local$parallel, ncores = local$ncores
        )
        de_ie <- spcov_params_val[["de"]] + spcov_params_val[["ie"]]
        c0 <- bq$c0
        s0 <- bq$s0 + de_ie * (1 / G - 1 / local$size_new)
      } else {
        # method_new = "basis": c0 (and x0) exact by parallelized (i.e., chunked) accumulation; s0
        # h(s) averaged exactly over every
        # block point, its outer average taken over the size_new nodes
        bq <- get_block_quantities(object, newdata,
          nodes = nodes,
          parallel = local$parallel, ncores = local$ncores
        )
        c0 <- bq$c0
        s0 <- bq$s0
      }
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
        keep <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, local$size))]
      } else if (local$method == "covariance") {
        # use abs() here for negative covariance types
        keep <- order(abs(as.numeric(c0)))[seq(from = n, to = max(1, n - local$size + 1))]
      }
      obdata <- obdata[keep, , drop = FALSE]
      c0 <- c0[keep]
      Xmat <- Xmat[keep, , drop = FALSE]
      y <- y[keep] # already offset-free (see above)
      cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(Sig[keep, keep, drop = FALSE])))
    }

    SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
    SqrtSigInv_y <- forwardsolve(cov_lowchol, y)
    residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
    SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)

    fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
    if (!is.null(offset_block)) {
      fit <- fit + offset_block
    }

    if (se.fit || interval == "prediction") {
      H <- x0 - crossprod(SqrtSigInv_c0, SqrtSigInv_X)
      vars <- as.numeric(s0 - crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% tcrossprod(cov_betahat, H))
      # a near-perfectly-predicted block can push the subtraction numerically below
      # zero; floor at 0 so sqrt() does not return NaN
      vars <- pmax(vars, 0)
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
    # apply the block's own (averaged) offset
    if (!is.null(offset_block)) {
      fit <- fit + offset_block
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

#' Accumulate the grid-dependent block-prediction covariance quantities
#'
#' A block prediction needs, over the \eqn{G} rows of the prediction grid,
#' \code{c0} (each observed site's average covariance with the block) and
#' \code{s0} (\eqn{Var(\bar{Y})} for the block; the average pairwise
#' covariance among the block points). Formed directly these are a dense
#' \eqn{G \times n_{obs}} and a dense \eqn{G \times G} matrix; here the grid
#' is passed in row-chunks and only column sums are kept, so memory is preserved.
#'
#' \code{s0} is written as a nested average,
#' \deqn{s0 = \frac{1}{G}\sum_i h(s_i), \qquad
#'       h(s) = \frac{1}{G}\sum_j Cov(s, s_j),}
#' and estimated by averaging \eqn{h}, computed exactly over all \eqn{G}
#' block points, at the rows named by \code{nodes}. With
#' \code{nodes = seq_len(G)} every block point is a node and \code{s0} is
#' exact (this replaces the old row-by-row \code{get_bk_var()}) used by \code{local$method_new = "basis"}.
#' \code{covmatrix()}'s \code{"pred.obs"} type omits the independent error
#' variance for a point against itself, so one \code{ie} is added back per
#' node (mirroring the old \code{get_each_bk_meancov()}), which also places
#' the block-variance diagonal at exactly the \eqn{1/G} weight that the exact
#' \code{s0} carries.
#'
#' @param object A fitted \code{splm} model object (its \code{obdata} is the
#'   observed data).
#' @param grid The prediction grid (coordinates already anisotropy-transformed
#'   and predictors already validated by \code{predict_block_splm()}).
#' @param nodes Integer row indices into \code{grid} giving the outer-average
#'   nodes for \code{s0}. \code{seq_len(NROW(grid))} gives the exact block
#'   variance.
#' @param chunk Number of \code{grid} rows accumulated per step.
#' @param parallel,ncores If \code{parallel}, the chunking loop is
#'   spread over \code{ncores} workers.
#'
#' @return A list with \code{c0} (length \code{n_obs}) and \code{s0} (scalar).
#'
#' @noRd
get_block_quantities <- function(object, grid, nodes, chunk = 1000L, parallel = FALSE, ncores = NULL) {
  G <- NROW(grid)
  ie <- object$coefficients$spcov[["ie"]]
  n_obs <- NROW(object$obdata)

  # object_nodes supplies grid[nodes, ] as the "observed" side so that
  # covmatrix(cov_type = "pred.obs") returns Cov(grid chunk, nodes), the
  # same obdata substitution covmatrix(cov_type = "pred.pred") uses
  object_nodes <- object
  object_nodes$obdata <- grid[nodes, , drop = FALSE]

  # keep each covariance block in memory
  chunks <- split(seq_len(G), ceiling(seq_len(G) / chunk))

  per_chunk <- function(rows) {
    grid_chunk <- grid[rows, , drop = FALSE]
    list(
      c0 = colSums(covmatrix(object, newdata = grid_chunk, cov_type = "pred.obs")),
      h = colSums(covmatrix(object_nodes, newdata = grid_chunk, cov_type = "pred.obs"))
    )
  }

  if (parallel) {
    ncores <- min(ncores, detectCores(), length(chunks))
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    parts <- parLapply(cl, chunks, per_chunk)
  } else {
    parts <- lapply(chunks, per_chunk)
  }

  c0 <- Reduce(`+`, lapply(parts, function(x) x$c0 / G))
  # + ie per node: the one self-covariance in each node's column that
  # "pred.obs" returned without the nugget, add it once to each row then average
  h <- (Reduce(`+`, lapply(parts, function(x) x$h / G)) + ie / G)

  list(c0 = as.numeric(c0), s0 = mean(h))
}
