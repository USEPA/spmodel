# shared building blocks for the predict.splm()/predict.spautor()/
# predict.spglm()/predict.spgautor() family (R/predict.R, R/predict_glm.R),
# plus conditional.R, decorrelate_newdata.R, and predict_block.R, which build
# prediction design matrices the same way.
#' Check that \code{newdata} contains the coordinate columns a point-referenced model needs
#'
#' \code{splm()}/\code{spglm()} point and block predictions all eventually
#' compute distances via \code{spdist_vectors()}, which treats a missing
#' coordinate column as a 0-length dimension rather than erroring immediately
#' -- the resulting mismatched or degenerate distance matrices only fail
#' several steps later, deep inside a Cholesky-related linear algebra call,
#' with an error that says nothing about coordinates or \code{newdata} (e.g.
#' "non-conformable arrays" or "invalid 'k' argument"; see
#' \code{\link{spdist_vectors}()}, which also guards against this
#' independently). Checking here, right after \code{newdata} is resolved
#' (sf-to-data-frame converted and 1D-zero-filled, if applicable) and before
#' any coordinate-dependent computation begins, gives a message that actually
#' explains what is wrong. \code{spautor()}/\code{spgautor()} have no
#' equivalent check -- their predictions use a precomputed neighbor structure
#' (\code{object$W}) rather than \code{newdata} coordinates.
#'
#' @param newdata The resolved \code{newdata}
#' @param xcoord,ycoord The x-coordinate/y-coordinate variable names
#'
#' @return Invisibly \code{NULL}; called for its error-checking side effect
#'
#' @noRd
check_newdata_coords <- function(newdata, xcoord, ycoord) {
  missing_coords <- setdiff(c(xcoord, ycoord), names(newdata))
  if (length(missing_coords) > 0) {
    stop(
      "newdata is missing the coordinate column(s) used to fit the model: ",
      paste0("\"", missing_coords, "\"", collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Build the prediction design matrix, guarding against a \code{poly()}/\code{model.matrix()} bug with one-row \code{newdata}
#'
#' When a formula term produces matrix-valued columns (e.g. \code{poly()},
#' \code{splines::bs()}/\code{ns()}) and \code{newdata} has exactly one row,
#' \code{stats::poly()}'s two-variable form (\code{poly(x, y, degree = n)})
#' misparses \code{y}'s single value as a positional \code{degree} argument
#' (see \code{stats::poly()}'s \code{dots_deg} heuristic), which either errors
#' outright or silently builds the wrong basis (I believe this is a bug with
#' base R that is hopefully fixed at some point). The single row is duplicated
#' to give \code{poly()} enough rows to build its basis correctly, and the
#' duplicate is dropped again once the model matrix is built. Shared by
#' \code{get_prediction_object_splm()}/\code{get_prediction_object_spautor()}
#' (R/get_prediction_object.R), their GLM counterparts
#' (R/get_prediction_object_glm.R), \code{conditional.splm()}/
#' \code{conditional.spglm()} (R/conditional.R), \code{decorrelate_newdata()}
#' (R/decorrelate_newdata.R), and \code{predict_block()} (R/predict_block.R).
#'
#' @param object A fitted model object
#' @param newdata The data requiring prediction
#'
#' @return A list with \code{newdata} (possibly reduced back to its original
#'   single row), \code{newdata_model}, and \code{offset} (from
#'   \code{model.offset()}, \code{NULL} if the formula has no \code{offset()}
#'   term)
#'
#' @noRd
get_newdata_model_matrix <- function(object, newdata) {
  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  # workaround: duplicate the single row so poly() has enough rows to build
  # its basis without erroring, build the model matrix, then keep only row 1
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
  list(newdata = newdata, newdata_model = newdata_model, offset = offset)
}

#' Resolve the extra random effect / partition factor lookup lists used by prediction
#'
#' \code{...} lets \code{loocv()} (which repeatedly calls \code{predict()} on
#' the same fitted object) pass a precomputed \code{extra_randcov_list}/
#' \code{extra_partition_list} through to avoid rebuilding these once per
#' left-out observation; ordinary \code{predict()} calls fall through to
#' building them fresh.
#'
#' @param object A fitted model object from \code{splm()}/\code{spglm()}
#' @param obdata The observed data
#' @param newdata The data requiring prediction
#' @param dotlist The calling predict method's \code{...}, as a named list
#'   (\code{list(...)})
#'
#' @return A list with \code{randcov_terms}, \code{reform_bar2}, and
#'   \code{partition_index_obdata}
#'
#' @noRd
get_extra_randcov_partition_lists <- function(object, obdata, newdata, dotlist) {
  dotlist_names <- names(dotlist)

  if ("extra_randcov_list" %in% dotlist_names && !is.null(dotlist[["extra_randcov_list"]])) {
    extra_randcov_list <- dotlist$extra_randcov_list
  } else {
    extra_randcov_list <- get_extra_randcov_list(object, obdata, newdata)
  }
  randcov_terms <- extra_randcov_list$randcov_terms

  if ("extra_partition_list" %in% dotlist_names && !is.null(dotlist[["extra_partition_list"]])) {
    extra_partition_list <- dotlist$extra_partition_list
  } else {
    extra_partition_list <- get_extra_partition_list(object, obdata, newdata)
  }
  reform_bar2 <- extra_partition_list$reform_bar2
  partition_index_obdata <- extra_partition_list$partition_index_obdata

  list(randcov_terms = randcov_terms, reform_bar2 = reform_bar2, partition_index_obdata = partition_index_obdata)
}

#' Precompute medium-sized-data prediction matrices, when it makes sense
#'
#' When using a distance/covariance local approximation, precomputing the
#' full observed-by-prediction distance, random effect, partition, and
#' covariance matrices once (vectorized across every prediction row) is
#' cheaper than having \code{get_pred_splm()}/\code{get_pred_spglm()}
#' reconstruct them from scratch for every row individually -- but only worth
#' it (a) when that matrix itself isn't too large to build and store in
#' memory, and (b) when a random effect or partition factor is present.
#' Confirmed empirically: without either, \code{get_pred_splm()}'s/
#' \code{get_pred_spglm()}'s per-row \code{spdist_vectors()}/\code{cov_vector()}
#' calls are already just as fast as this precompute (both are cheap,
#' vectorized operations with no meaningful per-call overhead) -- the real
#' cost this precompute amortizes is \code{randcov_vector()}'s/
#' \code{partition_vector()}'s underlying \code{Matrix::sparseMatrix()}
#' construction, which has real fixed per-call overhead in R's Matrix
#' package. So this is skipped entirely whenever there's no random
#' effect/partition factor to amortize it for, and its size is capped via
#' \code{local$byrow_threshold} (see \code{predict.spmodel()}) for users who
#' want to trade memory for speed differently than the default. Only used by
#' \code{predict.splm()}/\code{predict.spglm()} -- \code{predict.spautor()}/
#' \code{predict.spgautor()} have no analogous big-data approximation path.
#'
#' @param object A fitted model object from \code{splm()}/\code{spglm()}
#' @param newdata The data requiring prediction
#' @param obdata The observed data
#' @param xcoord,ycoord The x-coordinate/y-coordinate variable names
#' @param npred Number of prediction rows (\code{NROW(newdata)})
#' @param local_list A local prediction settings list (see
#'   \code{get_local_list_prediction()})
#' @param randcov_params_val A \code{randcov_params} object, or \code{NULL}
#' @param randcov_terms Random effect term context, from
#'   \code{get_extra_randcov_partition_lists()}
#' @param spcov_params_val A \code{spcov_params} object
#' @param reform_bar2,partition_index_obdata Partition factor context, from
#'   \code{get_extra_randcov_partition_lists()}
#'
#' @return A list with \code{dist_matrix_full}, \code{partition_vector_full},
#'   and \code{cov_vector_full} -- all \code{NULL} when the precompute is
#'   skipped
#'
#' @noRd
get_medium_precompute <- function(object, newdata, obdata, xcoord, ycoord, npred, local_list,
                                   randcov_params_val, randcov_terms, spcov_params_val,
                                   reform_bar2, partition_index_obdata) {
  if (local_list$method != "all" && (!is.null(object$random) || !is.null(object$partition_factor)) &&
    as.double(object$n) * npred < local_list$byrow_threshold) {
    dist_matrix_full <- spdist_vectors(newdata, obdata, xcoord, ycoord, object$dim_coords, sparse = FALSE)
    randcov_vector_full <- if (is.null(randcov_params_val)) {
      NULL
    } else {
      randcov_vector(randcov_params_val, obdata, newdata, randcov_terms)
    }
    partition_vector_full <- partition_vector(object$partition_factor, obdata, newdata, reform_bar2, partition_index_obdata)
    cov_vector_full <- cov_vector(spcov_params_val, dist_matrix_full, randcov_vector_full, partition_vector_full)
  } else {
    dist_matrix_full <- NULL
    partition_vector_full <- NULL
    cov_vector_full <- NULL
  }
  list(dist_matrix_full = dist_matrix_full, partition_vector_full = partition_vector_full, cov_vector_full = cov_vector_full)
}

#' Dispatch per-row prediction, in parallel or sequentially
#'
#' @param fn The per-row prediction function (\code{get_pred_splm()},
#'   \code{get_pred_spautor()}, \code{get_pred_spglm()}, or
#'   \code{get_pred_spgautor()})
#' @param data_list A list with one element per prediction row (either
#'   \code{newdata} split by row, or, for the areal methods, a cluster list of
#'   per-row \code{x0}/\code{c0}/\code{s0} quantities)
#' @param prediction_object The prediction context object, spliced with
#'   \code{pred_row_context} and passed to \code{fn} unchanged for every row
#' @param local_list A local prediction settings list (see
#'   \code{get_local_list_prediction()}) -- \code{local_list$parallel} selects
#'   \code{parallel::parLapply()} over \code{local_list$ncores} workers versus
#'   plain \code{lapply()}
#'
#' @return A list with one element per prediction row, as returned by \code{fn}
#'
#' @noRd
run_pred_dispatch <- function(fn, data_list, prediction_object, local_list) {
  if (local_list$parallel) {
    cl <- parallel::makeCluster(local_list$ncores)
    pred_list <- parallel::parLapply(cl, data_list, fn, prediction_object = prediction_object)
    cl <- parallel::stopCluster(cl)
  } else {
    pred_list <- lapply(data_list, fn, prediction_object = prediction_object)
  }
  pred_list
}

#' Apply row naming and pick the return shape for an \code{interval = "none"} prediction
#'
#' The last step of every \code{interval == "none"} branch across
#' \code{predict.splm()}/\code{predict.spautor()}/\code{predict.spglm()}/
#' \code{predict.spgautor()} is identical once \code{fit}/\code{se} are on
#' their final scale (post-offset, post-\code{invlink()}/delta-method
#' adjustment where applicable) -- name the rows if requested and return
#' either \code{fit} alone or \code{list(fit, se.fit)}, depending on whether
#' standard errors were requested at all.
#'
#' @param fit The point predictions, already on their final scale
#' @param se The standard errors, already on their final scale, or
#'   \code{NULL} when standard errors were not requested (\code{se.fit ==
#'   FALSE}) -- this \code{NULL}-ness is itself the signal for which return
#'   shape to use, since \code{se} (unlike in
#'   \code{finalize_interval_bounds()}) is never computed at all when not
#'   requested
#' @param add_newdata_rows Whether to name the returned values using
#'   \code{missing_index} -- always \code{TRUE} for \code{predict.spautor()}/
#'   \code{predict.spgautor()}, since areal prediction has no separate
#'   arbitrary-\code{newdata} concept to guard against
#' @param missing_index Row labels to apply (\code{object$missing_index})
#'
#' @return \code{list(fit = fit, se.fit = se)} if \code{se} was supplied,
#'   otherwise \code{fit} alone
#'
#' @noRd
finalize_interval_none <- function(fit, se, add_newdata_rows, missing_index) {
  if (!is.null(se)) {
    if (add_newdata_rows) {
      names(fit) <- missing_index
      names(se) <- missing_index
    }
    list(fit = fit, se.fit = se)
  } else {
    if (add_newdata_rows) {
      names(fit) <- missing_index
    }
    fit
  }
}

#' Assemble the fit/lwr/upr matrix and pick the return shape for an interval prediction
#'
#' The last step of every \code{interval == "prediction"}/\code{"confidence"}
#' branch across all four \code{predict.*()} methods is identical once
#' \code{fit}/\code{lwr}/\code{upr}/\code{se} are on their final scale
#' (post-offset, post-\code{invlink()}/delta-method adjustment where
#' applicable) -- bind them into the returned matrix, name the rows if
#' requested, and return either the matrix alone or
#' \code{list(fit, se.fit)}. Unlike \code{finalize_interval_none()},
#' \code{se} is always a real vector here (needed upstream to build
#' \code{lwr}/\code{upr} regardless of whether standard errors were actually
#' requested), so \code{se.fit} must be passed explicitly rather than
#' inferred from \code{se}'s \code{NULL}-ness.
#'
#' @param fit,lwr,upr The point predictions and interval bounds, already on
#'   their final scale
#' @param se The standard errors, already on their final scale
#' @param se.fit Whether standard errors were requested
#' @param add_newdata_rows Whether to name the returned values using
#'   \code{missing_index} -- always \code{TRUE} for \code{predict.spautor()}/
#'   \code{predict.spgautor()}, since areal prediction has no separate
#'   arbitrary-\code{newdata} concept to guard against
#' @param missing_index Row labels to apply (\code{object$missing_index})
#'
#' @return \code{list(fit = <n x 3 matrix>, se.fit = se)} if \code{se.fit},
#'   otherwise the \code{<n x 3 matrix>} alone
#'
#' @noRd
finalize_interval_bounds <- function(fit, lwr, upr, se, se.fit, add_newdata_rows, missing_index) {
  fit <- cbind(fit, lwr, upr)
  row.names(fit) <- seq_len(NROW(fit))
  if (se.fit) {
    if (add_newdata_rows) {
      row.names(fit) <- missing_index
      names(se) <- missing_index
    }
    list(fit = fit, se.fit = se)
  } else {
    if (add_newdata_rows) {
      row.names(fit) <- missing_index
    }
    fit
  }
}

#' Build the local distance/covariance vector for a single \code{splm()}/\code{spglm()} prediction row
#'
#' Shared body of \code{get_pred_splm()}/\code{get_pred_spglm()} up through
#' the distance/covariance vector needed for nearest-neighbor/nearest-
#' covariance subsetting -- handles the medium-mode data reuse (slicing this
#' row out of the full observed-by-prediction matrices \code{predict.splm()}/
#' \code{predict.spglm()} already built, when available), the partition-index
#' subsetting, and the dense-mode from-scratch recompute. Everything after
#' this point (the actual nearest-neighbor/nearest-covariance subsetting, the
#' local Cholesky build, and whitening) differs by family -- Gaussian works on
#' \code{y}, recomputed from \code{model.frame()} when locally subsetting; GLM
#' works on \code{w}/\code{size}, which must be subsetted directly since they
#' aren't recomputed from a formula -- so it stays inline in each caller.
#'
#' @param newdata_list A single element of \code{newdata} split by row
#'   (\code{newdata_list$row}), with \code{newdata_list$row_index} present
#'   only in medium mode
#' @param obdata The observed data (possibly already partition-subsetted by
#'   the caller)
#' @param xcoord,ycoord The x-coordinate/y-coordinate variable names
#' @param dim_coords Number of coordinate dimensions
#' @param spcov_params_val A \code{spcov_params} object
#' @param randcov_params_val A \code{randcov_params} object, or \code{NULL}
#' @param randcov_terms Random effect term context
#' @param partition_factor,reform_bar2,partition_index_obdata Partition
#'   factor context
#' @param random The random effect formula, or \code{NULL}
#' @param local A local prediction settings list (\code{local$method} selects
#'   \code{"all"}/\code{"distance"}/\code{"covariance"})
#' @param dist_matrix_full,partition_vector_full,cov_vector_full The
#'   medium-mode precomputed matrices from \code{get_medium_precompute()}, or
#'   all \code{NULL} when medium mode is not in use
#'
#' @return A list with \code{obdata}, \code{randcov_terms} (both possibly
#'   partition-subsetted), \code{dist_vector}, and \code{cov_vector_val}
#'
#' @noRd
get_pred_local_setup <- function(newdata_list, obdata, xcoord, ycoord, dim_coords,
                                  spcov_params_val, randcov_params_val, randcov_terms,
                                  partition_factor, reform_bar2, partition_index_obdata,
                                  random, local, dist_matrix_full, partition_vector_full,
                                  cov_vector_full) {
  # medium-sized data: the caller already built the full observed-by-
  # prediction distance/covariance matrices vectorized across every
  # prediction row, so slice this row out of them instead of reconstructing
  # from scratch (identical values, since matrix entries don't depend on
  # whether they were computed one row at a time or all at once)
  medium_mode <- !is.null(cov_vector_full)

  if (medium_mode) {
    row_index <- newdata_list$row_index
    partition_vector <- if (is.null(partition_vector_full)) NULL else partition_vector_full[row_index, , drop = FALSE]
    dist_vector <- dist_matrix_full[row_index, , drop = FALSE]
    cov_vector_val <- cov_vector_full[row_index, , drop = FALSE]
  } else {
    # storing partition vector
    partition_vector <- partition_vector(partition_factor,
      data = obdata,
      newdata = newdata_list$row, reform_bar2 = reform_bar2,
      partition_index_data = partition_index_obdata
    )
  }

  # subsetting partition vector (efficient but causes problems later with
  # random effect subsetting)
  # this must happen BEFORE the distance/covariance neighbor search further
  # below: observations outside the new observation's partition group have
  # zero covariance with it by construction, so if the local$size nearest/
  # most-correlated neighbors were instead selected from the full (un-
  # partitioned) obdata, method = "distance" -- which ranks by raw spatial
  # distance and has no notion of partition membership -- could silently pick
  # neighbors from another partition, and method = "covariance" would waste
  # neighbor slots on observations whose covariance is always exactly zero
  if (!is.null(partition_vector) && local$method %in% c("distance", "covariance") &&
    (is.null(random) || !labels(terms(partition_factor)) %in% labels(terms(random)))) {
    partition_index <- as.vector(partition_vector) == 1
    randcov_terms <- lapply(randcov_terms, function(term) {
      term$group_label <- term$group_label[partition_index]
      # level_index_map indexes into group_label, so it must be rebuilt against
      # the subsetted vector rather than reused from the full obdata
      term$level_index_map <- split(seq_along(term$group_label), term$group_label)
      if (!is.null(term$slope_val)) {
        term$slope_val <- term$slope_val[partition_index]
      }
      term
    })
    obdata <- obdata[partition_index, , drop = FALSE]
    partition_vector <- Matrix(1, nrow = 1, ncol = NROW(obdata))
    # medium-mode's dist_vector/cov_vector_val were sliced from the full
    # observed-by-prediction matrices above, so they must be restricted to
    # the same partition_index columns to stay aligned with the now-
    # subsetted obdata
    if (medium_mode) {
      dist_vector <- dist_vector[, partition_index, drop = FALSE]
      cov_vector_val <- cov_vector_val[, partition_index, drop = FALSE]
    }
  }

  # dense (non-medium) path: recompute the distance/covariance vectors from
  # scratch against the (possibly partition-subsetted) obdata, since medium
  # mode's precomputed full matrices aren't available here
  if (!medium_mode) {
    # dense: a distance vector/matrix has essentially no exact zeros, so sparse
    # storage buys nothing but S4 dispatch overhead on every prediction row
    dist_vector <- spdist_vectors(newdata_list$row, obdata, xcoord, ycoord, dim_coords, sparse = FALSE)

    # making random vector if necessary
    if (!is.null(randcov_params_val)) {
      randcov_vector_val <- randcov_vector(randcov_params_val, obdata, newdata_list$row, randcov_terms)
    } else {
      randcov_vector_val <- NULL
    }

    # making the covariance vector
    cov_vector_val <- cov_vector(spcov_params_val, dist_vector, randcov_vector_val, partition_vector)
  }

  list(obdata = obdata, randcov_terms = randcov_terms, dist_vector = dist_vector, cov_vector_val = cov_vector_val)
}

#' Compute the fit/var prediction for a single areal (autoregressive) cluster
#'
#' Shared body of \code{get_pred_spautor()}/\code{get_pred_spgautor()} --
#' identical math for both families; only \code{residuals_pearson}'s
#' *meaning* differs (observed-scale Pearson residuals for \code{splm}'s
#' Gaussian response vs. link-scale residuals for \code{spglm}'s GLM latent
#' variable), not its role in the computation. Unlike \code{get_pred_splm()}/
#' \code{get_pred_spglm()}, there is no medium-mode/local-neighbor-search
#' distinction here -- areal models use one fixed joint covariance matrix
#' over all sites, already sliced down to this cluster by the caller.
#'
#' @param cluster_list A list with elements \code{x0} (the design matrix row
#'   for the new observation), \code{c0} (the covariance vector between the
#'   new observation and the observed data), and \code{s0} (the marginal
#'   variance of the new observation)
#' @param cov_matrix_lowchol The lower Cholesky factor of the observed-data
#'   covariance matrix
#' @param betahat,cov_betahat The estimated coefficients and their covariance
#' @param residuals_pearson The observed-data Pearson (or link-scale, for GLM)
#'   residuals
#' @param SqrtSigInv_X \code{cov_matrix_lowchol^{-1} \%*\% Xmat}
#' @param se.fit,interval,type As in \code{predict.spautor()}/\code{predict.spgautor()}
#' @param Xmat The observed-data design matrix
#'
#' @return A list with element \code{fit} (and, if \code{se.fit} or
#'   \code{interval == "prediction"}, element \code{var}) for the new observation
#'
#' @noRd
get_areal_pred <- function(cluster_list, cov_matrix_lowchol, betahat, residuals_pearson,
                                cov_betahat, SqrtSigInv_X, se.fit, interval, type, Xmat) {
  x0 <- cluster_list$x0
  c0 <- cluster_list$c0
  s0 <- cluster_list$s0
  SqrtSigInv_c0 <- forwardsolve(cov_matrix_lowchol, c0)
  if (type == "weight") {
    Xt_SigInv <- t(backsolve(t(cov_matrix_lowchol), SqrtSigInv_X))
    betahat_wt <- cov_betahat %*% Xt_SigInv
    residuals_weight <- -1 * Xmat %*% betahat_wt # this is recomputed over and over when using all data consider making more efficient
    diag(residuals_weight) <- diag(residuals_weight) + 1
    fit <- x0 %*% betahat_wt + Matrix::crossprod(SqrtSigInv_c0, forwardsolve(cov_matrix_lowchol, residuals_weight))
  } else {
    # universal kriging BLUP: trend plus a covariance-weighted combination of
    # the observed residuals
    fit <- as.numeric(x0 %*% betahat + crossprod(SqrtSigInv_c0, residuals_pearson))
  }

  if (se.fit || interval == "prediction") {
    # kriging variance decomposition: marginal variance s0, minus variance
    # explained by the observed data, plus the betahat-estimation-uncertainty
    # inflation term
    H <- x0 - crossprod(SqrtSigInv_c0, SqrtSigInv_X)
    var <- as.numeric(s0 - crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% tcrossprod(cov_betahat, H))
    pred_list <- list(fit = fit, var = var)
  } else {
    pred_list <- list(fit = fit)
  }
  pred_list
}
