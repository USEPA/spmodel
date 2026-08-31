#' Conditionally simulate from a model
#'
#' @description Conditionally simulate prediction data from a model object.
#'
#' @param object A fitted model object.
#' @param newdata A data frame or \code{sf} object in which to
#'   look for variables with which to predict. If a data frame, \code{newdata}
#'   must contain all variables used by \code{formula(object)} and all variables
#'   representing coordinates. If an \code{sf} object, \code{newdata} must contain
#'   all variables used by \code{formula(object)} and coordinates are obtained
#'   from the geometry of \code{newdata}. If omitted, missing data from the
#'   fitted model object are used.
#' @param output  The output type, which can be any subset of
#'   \code{c("newdata", "beta", "object")}. If \code{simulate_covparams = TRUE},
#'   the output type can also be any subset of \code{c("cov", "spcov", "randcov")}.
#'    The default is \code{"newdata"}. See Details for more.
#' @param type For \code{spglm()} model objects, the scale of the conditional
#'   simulations for \code{newdata}.
#'   When \code{type = "link"}, the predicted means on the
#'   link scale are returned. When \code{type = "response"}, the predicted means
#'   on the response scale are returned. When \code{type = "new"}, a new observation
#'   is simulated from the appropriate response distribution with mean equal to
#'   the mean on the response scale and dispersion equal to the dispersion parameter
#'   from \code{object}. The default is \code{"link"}.
#' @param samples  The number of conditional simulations. The default is
#'   \code{1,000}.
#' @param local An optional logical or list controlling the big data approximation.
#'   If omitted, \code{local} is set
#'   to \code{TRUE} or \code{FALSE} based on the whether the observed
#'   or prediction sample size (the number of
#'   non-missing observations in \code{data} or \code{newdata}) exceeds 5,000,
#'   \code{local} is set to \code{TRUE}. Otherwise it is set to \code{FALSE}.
#'   If \code{FALSE}, no big data approximation is implemented.
#'   If a list is provided, \code{local$approximation} selects which big data
#'   approximation is used and can take on the values
#'   \code{"low-rank"} or \code{"vecchia"}:
#'   \itemize{
#'     \item \code{"low-rank"}: a base sample is drawn from the observed
#'       data, \code{newdata} is split into blocks, and each block is
#'       simulated conditional on the base sample alone  (blocks are
#'       assumed conditionally independent of one another given the base
#'       sample). The base-sample settings (\code{method_base}/\code{size_base}/
#'       \code{reorder_base}) and the \code{newdata}-blocking settings
#'       (\code{method_new}/\code{size_new}/\code{reorder_new}/\code{kmeans_new})
#'       are separate from one another.
#'       \itemize{
#'         \item \code{method_base}: Whether the observed data conditioned on is
#'           restricted to a base sample. If \code{method_base = "all"}, no
#'           big data approximation is applied to the observed data (all of it is
#'           conditioned on). If \code{method_base = "base"}, the observed data
#'           is subset to \code{size_base} observations (ordered via
#'           \code{reorder_base}) to form the base sample. The default is
#'           \code{"base"}.
#'         \item \code{reorder_base}: The data reordering approached to reorder
#'           the observed data prior to subsetting to obtain a base sample.
#'           If \code{reorder = "none"}, no reordering
#'           is applied to the observed data. If \code{reorder = "random"}, the observed data order is
#'           randomly reshuffled. If \code{reorder = "grts"}, the observed data order is
#'           randomly generated using the GRTS algorithm for spatially balanced
#'           sampling via \code{spsurvey::grts()}. The default is \code{"grts"}.
#'         \item \code{size_base}: The number of observed data observations used for the base sample.
#'           The default is 5,000. See Details for more.
#'         \item \code{method_new}: Whether \code{newdata} is split into blocks
#'           for simulation. If \code{method_new = "all"}, no big data
#'           approximation is applied to \code{newdata} (it is simulated all at
#'           once). If \code{method_new = "base"}, \code{newdata} is split into
#'           blocks of (approximately) \code{size_new} observations each (ordered
#'           via \code{reorder_new}, and optionally grouped via \code{kmeans_new}),
#'           with each block simulated conditional on the base sample
#'           independently of every other block. The default is \code{"base"}.
#'         \item \code{reorder_new}: The data reordering approach used to reorder
#'           \code{newdata} before splitting it into blocks. If \code{reorder = "none"}, no reordering
#'           is applied to \code{newdata}. If \code{reorder = "random"}, \code{newdata} is
#'           randomly reshuffled. The default is \code{"random"}.
#'         \item \code{kmeans_new}: Whether \code{newdata} observations
#'           should be assigned to blocks based on k-means clustering
#'           on the coordinates, with clusters of size approximately equal to
#'           \code{size_new}. The default is \code{FALSE} when \code{reorder_new = "none"}
#'           and \code{TRUE} otherwise.
#'         \item \code{size_new}: The (approximate) number of observations used
#'           for each block. The default is 1,000. See Details for more.
#'         \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'           parallel package is automatically used. The default is \code{FALSE}.
#'         \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'           parallelize over. The default is the number of available cores on your machine.
#'       }
#'       If \code{local$approximation} is \code{"low-rank"} (either explicitly or via
#'       \code{local = TRUE}), defaults for the remaining \code{"low-rank"}
#'       settings are chosen such that \code{local} is transformed into
#'       \code{list(approximation = "low-rank", method_base = "base", size_base = 5000,
#'       reorder_base = "grts", method_new = "base", size_new = 1000,
#'       reorder_new = "random", kmeans_new = TRUE, parallel = FALSE)}.
#'     \item \code{"vecchia"}: every \code{newdata} location is simulated one
#'       at a time (in some order over \code{newdata}), each conditional on
#'       \strong{all} observed data plus every already-simulated
#'       \code{newdata} location (not a single shared base sample).
#'       \code{newdata} locations are never assumed conditionally independent
#'       of one another. This is exact (matches \code{local = FALSE}) when
#'       \code{method = "all"}; \code{method = "distance"}/\code{"covariance"}
#'       truncate the conditioning set to a fixed number of neighbors
#'       sorted by distance or covariance with the new observation. No parallelization
#'       exists because the algorithm is inherently sequential, as each new observation
#'       depends on previous ones.
#'       \itemize{
#'         \item \code{method}: The neighbor-selection rule used to build each
#'           location's conditioning set once it exceeds \code{size} candidates
#'           (all observed data plus every already-simulated \code{newdata}
#'           location). Values are \code{"all"}, \code{"distance"}
#'           (the \code{size} nearest candidates), or \code{"covariance"} (the
#'           \code{size} candidates with the highest covariance, in absolute
#'           value, with the location being simulated). Same convention as
#'           \code{predict()}'s own \code{local$method}. The default is
#'           \code{"covariance"}. \code{method = "all"} is very computationally
#'           intensive and \code{local = FALSE} should almost always be used instead.
#'           (\code{method = "all"} primarily exists for numerical verification).
#'         \item \code{size}: The number of neighbors used when \code{method}
#'           is \code{"distance"} or \code{"covariance"}. The default is 30.
#'         \item \code{ordering}: The order \code{newdata} locations are
#'           simulated in -- \code{"maxmin"}, \code{"middleout"},
#'           \code{"outsidein"}, \code{"coordinate"}, \code{"grts"},
#'           \code{"random"}, or \code{"none"} (same options as
#'           \code{decorrelate()}'s \code{ordering} argument). The default is
#'           \code{"maxmin"}.
#'       }
#'   }
#'       When \code{local = TRUE}, \code{local} is transformed into
#'       \code{list(approximation = "low-rank", method_base = "base", size_base = 5000,
#'       reorder_base = "grts", method_new = "base", size_new = 1000,
#'       reorder_new = "random", kmeans_new = TRUE, parallel = FALSE)}.
#'       When \code{local} is a list, at least one list element must be provided to
#'       initialize default arguments for the other list elements. See Details for more.
#' @param simulate_covparams For \code{splm()} model objects, whether to also
#'   simulate new covariance parameters for each sample. \code{simulate_covparams}
#'   requires \code{object$vcov$cov} to be specified during model fitting by
#'   selecting \code{ddf = "satterthwaite"}. \code{simulate_covparams = TRUE}
#'   should generally not be used for sample sizes greater than 500
#'   given its computational inefficiencies. The default is
#'   \code{FALSE}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#' @param newdata_size The \code{size} value for each observation in \code{newdata}
#'   used when predicting for the binomial family, with a default value of 1.
#'
#' @details
#'
#'   If \code{"newdata"} is in \code{output},
#'   conditional simulations are returned for each row of \code{newdata}.
#'   If \code{"beta"} is in \code{output},
#'   conditional simulations are returned for each fixed effect
#'   (i.e., element of \code{coef(object)}. If \code{"object"} is in \code{output},
#'   the observed data from \code{object} is returned once for each row of
#'   \code{newdata}. For example, \code{c("newdata", "beta")} returns
#'   the conditional simulations both for \code{newdata} and for the
#'   fixed effects. If \code{"cov"}/\code{"spcov"}/\code{"randcov"} is in
#'   \code{output} (only available when \code{simulate_covparams = TRUE}),
#'   the simulated covariance parameter draws themselves are returned.
#'
#'   \code{local} Details: When \code{local$approximation} is \code{"low-rank"}, the big
#'   data approximation works by assigning \code{size_base}
#'   observations to a base sample and then simulating data for the base sample.
#'   The remaining observations are assigned to blocks. For each block, data
#'   are simulated from the conditional distribution given the base sample.
#'   Observations from the same block share conditional covariance while
#'   observations from distinct blocks are assumed conditionally independent
#'   (given the base sample). Parallelization generally further speeds up
#'   computations. When \code{local$approximation} is \code{"vecchia"}, no such
#'   independence assumption is made -- see the \code{local} argument above
#'   for details. For \code{spglm()} model objects, both \code{local$approximation}s
#'   propagate the latent process's own estimation uncertainty
#'   (\code{var_adj}) analytically rather than by simulation; for
#'   \code{"vecchia"} this requires factorizing a dense matrix over all
#'   observed data one time, since this particular source of uncertainty is
#'   not spatially local and so cannot be shrunk by neighbor truncation the
#'   way the rest of the simulation is -- see the \code{local} argument above.
#'
#' @return If \code{output = "newdata"}, an a x b matrix of conditional simulations
#'   for each row in \code{newdata}, where a is the
#'   number of rows in \code{newdata} and b is the number of samples.
#'   If \code{output = "beta"}, an p x b matrix of conditional simulations for each
#'   element in \code{coef(object)}, where p is the
#'   number of fixed effects and b is the number of samples.
#'   If \code{output = "object"}, an n x b matrix of observed data values, where n is the
#'   number of rows in \code{data} and b is the number of samples.
#'   If \code{output = "cov"}/\code{"spcov"}/\code{"randcov"}
#'   (\code{simulate_covparams = TRUE} only), a (covariance parameter) x b
#'   matrix of the of conditoinal simulations for each covariance parameter, where b is the
#'   number of samples. \code{"cov"} returns every covariance parameter,
#'   while \code{"spcov"}/\code{"randcov"} return just the spatial/random-effect
#'   covariance parameters, respectively.
#'
#'   If \code{output} has more than one element, a list is returned with the
#'   respetive elements named according to the relevant \code{output}. For example
#'   \code{output = c("newdata", "beta")} returns a list with elements
#'   \code{"newdata"} and \code{"beta"}, each containing the relevant output
#'   for \code{output = "newdata"} and \code{output = "beta"}, respectively.
#'
#' @export
#'
#' @examples
#' set.seed(0)
#' spmod <- splm(sulfate ~ 1, data = sulfate, spcov_type = "exponential")
#' cond <- conditional(spmod, newdata = sulfate_preds)
#' predict(spmod, sulfate_preds[20, ], se.fit = TRUE)
#' c("fit_cond" = mean(cond[20, ]), "se.fit_cond" = sd(cond[20, ]))
#' hist(cond[20, ])
conditional <- function(object, ...) {
  UseMethod("conditional", object)
}

#' @rdname conditional
#' @method conditional splm
#' @export
conditional.splm <- function(object, newdata, output = "newdata", samples = 1000, local, simulate_covparams = FALSE, ...) {

  if (missing(local)) {
    local <- NULL
  }
  if (!is.logical(simulate_covparams) || length(simulate_covparams) != 1 || is.na(simulate_covparams)) {
    stop("simulate_covparams must be TRUE or FALSE.", call. = FALSE)
  }

  if ("all" %in% output) {
    output <- c("newdata", "beta", "object")
  }
  if (any(!output %in% c("newdata", "beta", "object", "cov", "spcov", "randcov"))) {
    stop("output must be \"newdata\", \"beta\", \"object\", \"cov\", \"spcov\", \"randcov\", or \"all\".", call. = FALSE)
  }

  # error if newdata missing from arguments and object
  if (missing(newdata)) {
    if (is.null(object$newdata)) {
      stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
    } else {
      newdata <- object$newdata
    }
  }

  # local_list/simulate_covparams/samples must all be specified before samples is
  # first used below (the output = "object" shortcut)
  # get_local_list_conditional() only needs object/newdata
  local_list <- get_local_list_conditional(local, object, newdata)

  # simulate_covparams cannot reuse a single shared covariance factorization
  # across samples the way the rest of conditional() does, so it is only
  # supported when no big-data approximation is actually in effect. Once
  # local_list$method_base/method_new are both "all" (guaranteed here), the
  # simulate_covparams path below always operates on the full observed data
  # and all of newdata as one block, with no base-subsampling/blocking to
  # account for
  local_active <- local_list$approximation == "vecchia" ||
    (local_list$approximation == "low-rank" && (local_list$method_base != "all" || local_list$method_new != "all"))
  if (isTRUE(simulate_covparams) && local_active) {
    simulate_covparams <- FALSE
    message("simulate_covparams = TRUE is not used when a big-data approximation (local) is specified; setting simulate_covparams = FALSE.")
  }

  if (!isTRUE(simulate_covparams) && any(c("cov", "spcov", "randcov") %in% output)) {
    stop("output can only include \"cov\", \"spcov\", or \"randcov\" when simulate_covparams = TRUE.", call. = FALSE)
  }

  if (isTRUE(simulate_covparams)) {
    if (object$n > 500) {
      warning("simulate_covparams = TRUE may result in exceedingly long computational times when observed data sample sizes greater than 500.", call. = FALSE)
    }
    vcov_theta <- get_vcov_theta_for_conditional(object)
  }

  y <- model.response(model.frame(object))
  # output = "object" just replicates the observed y across simulation columns
  # (for easy row-binding with the newdata/beta draws below)
  base_val_y <- matrix(rep(y, times = samples), ncol = samples)
  if (length(output) == 1 && output == "object") {
    return(base_val_y)
  }
  # handle offset
  offset_obdata <- model.offset(model.frame(object))
  if (!is.null(offset_obdata)) {
    y <- y - offset_obdata
  }

  betahat <- coef(object)
  X <- model.matrix(object)

  # now simulate beta and add
  newdata_model_list <- get_newdata_model_matrix(object, newdata)
  newdata <- newdata_model_list$newdata
  newdata_model <- newdata_model_list$newdata_model
  offset <- newdata_model_list$offset
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  # keep only the newdata_model columns that also appear in the fitted
  # design matrix (factor levels present in data but absent from newdata can
  # otherwise leave newdata_model with a different column layout than X)
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  if (isTRUE(simulate_covparams)) {
    # covariance parameters resimulated per draw -- see get_conditional_covparams()
    # for why this cannot share the vectorized-across-samples approach below
    covparams_out <- get_conditional_covparams(object, newdata, y, X, betahat, vcov_theta, samples)
    new_betahat <- covparams_out$new_betahat
    new_val <- covparams_out$new_val
    new_cov <- covparams_out$new_cov
    new_spcov <- covparams_out$new_spcov
    new_randcov <- covparams_out$new_randcov
  } else {
    new_cov <- NULL
    new_spcov <- NULL
    new_randcov <- NULL
    # Composition sampling strategy for p(y0 | y): rather than drawing y0
    # directly from its fitted predictive
    # distribution, first draw beta from its asymptotic sampling distribution
    # N(betahat, vcov(object)) via a Cholesky factor, then (below) draw the
    # spatial residual field conditional on the observed residuals implied by
    # each drawn beta, and finally add the drawn beta's trend back in. This
    # propagates fixed effect uncertainty into the conditional draws instead of
    # conditioning on betahat alone.
    cov_betahat_lowchol <- t(chol(vcov(object)))
    new_betahat <- vapply(seq_len(samples), function(x) as.numeric(cov_betahat_lowchol %*% rnorm(length(betahat))), numeric(length(betahat)))
    # beta0 force to matrix
    if (!is.matrix(new_betahat)) {
      new_betahat <- matrix(new_betahat, nrow = 1)
    }
    new_betahat <- sweep(new_betahat, 1, betahat, "+")
    new_fitted <- X %*% new_betahat
    # residualize the observed y against each simulated beta -- these
    # (mean-zero, per-draw) residuals are what the spatial field is actually
    # conditioned on below; the simulated trend is added back in at the end
    new_resid <- sweep(-1 * new_fitted, 1, y, "+")

    if (local_list$approximation == "vecchia") {
      # vecchia: every newdata location is simulated sequentially, conditional
      # on ALL observed data (not subsampled, see get_conditional_vecchia())
      # plus every earlier-simulated newdata location, so there is no base/block
      # splitting step here at all
      new_val <- get_conditional_vecchia(object, newdata, new_resid, local_list, samples)
    } else {
      # low-rank, part 1: restrict the "observed" data conditioned on to a
      # spatially-representative base sample instead of all of data, so the
      # base covariance matrix factorized below stays a manageable size
      if (local_list$method_base != "all") {
        object$obdata <- object$obdata[local_list$index$base, , drop = FALSE]
        y_base <- y[local_list$index$base]
        X_base <- X[local_list$index$base, , drop = FALSE]
      } else {
        y_base <- y
        X_base <- X
      }
      # low-rank, part 2: split newdata into blocks so each block's
      # observed-by-prediction covariance is computed and factorized
      # separately (blocks are treated as conditionally independent given the
      # base sample); each block is a list element handled by
      # get_conditional_new_from_base_adjust() below, in parallel if requested
      if (local_list$method_new != "all") {
        x0 <- lapply(local_list$index$new, function(x) newdata_model[x, , drop = FALSE])
        newdata_split <- lapply(local_list$index$new, function(x) newdata[x, , drop = FALSE])
      } else {
        x0 <- list(newdata_model)
        newdata_split <- list(newdata)
      }
      newdata_list <- mapply(x = x0, y = newdata_split, FUN = function(x, y) list(x0 = x, newdata = y), SIMPLIFY = FALSE)
      spcov_val <- coef(object, type = "spcov")
      # pure nugget (independent error, no spatial dependence or random effects):
      # the base covariance matrix is diagonal, so a plain sqrt() gives its
      # (lower triangular) Cholesky factor without paying for a full chol()
      if (spcov_val[["de"]] == 0 && is.null(coef(object, type = "randcov"))) {
        # object$n is the full sample size not the base sample size
        cov_lowchol_base <- Matrix::Diagonal(n = NROW(object$obdata), x = sqrt(spcov_val[["ie"]]))
      } else {
        cov_lowchol_base <- t(chol(covmatrix(object)))
      }
      # solved once and shared across every block instead of re-solving a
      # full n_base x samples residual matrix inside each block's call (see
      # get_conditional_new_from_base_adjust())
      SqrtSigInv_y <- forwardsolve(cov_lowchol_base, cbind(y_base))
      SqrtSigInv_X <- forwardsolve(cov_lowchol_base, X_base)
      if (local_list$parallel) {
        cl <- parallel::makeCluster(local_list$ncores)
        new_val <- parLapply(cl, newdata_list, get_conditional_new_from_base_adjust, object, SqrtSigInv_y, SqrtSigInv_X, new_betahat, cov_lowchol_base, samples)
        cl <- parallel::stopCluster(cl)
      } else {
        new_val <- lapply(newdata_list, get_conditional_new_from_base_adjust, object, SqrtSigInv_y, SqrtSigInv_X, new_betahat, cov_lowchol_base, samples)
      }


      new_val <- do.call("rbind", new_val)
      # blocks were processed independently (and possibly reordered upstream by
      # GRTS/k-means grouping in get_local_list_conditional()), so restore the
      # original newdata row order before returning
      if (local_list$method_new != "all") {
        index_new <- do.call("c", local_list$index$new)
        new_val <- new_val[order(index_new), , drop = FALSE]
      }
    }
  }

  # add the simulated fixed effect trend (X0 %*% beta_b) back onto the
  # (mean-zero) conditional spatial residual draws -- completes the
  # composition sampling described previously
  new_val <- newdata_model %*% new_betahat + new_val

  if (!is.null(offset)) {
    new_val <- sweep(new_val, 1, offset, "+")
  }

  val <- list(newdata = new_val, beta = new_betahat, object = base_val_y, cov = new_cov, spcov = new_spcov, randcov = new_randcov)
  if (length(output) == 1) {
    return(val[[output]])
  } else {
    return(val[output])
  }
}

#' @rdname conditional
#' @method conditional spglm
#' @export
conditional.spglm <- function(object, newdata, output = "newdata", type = c("link", "response", "new"), samples = 1000, local, newdata_size, ...) {

  if (missing(local)) {
    local <- NULL
  }

  if ("all" %in% output) {
    output <- c("newdata", "beta", "object")
  }
  if (any(!output %in% c("newdata", "beta", "object"))) {
    stop("output must be \"newdata\", \"beta\", \"object\", or \"all\".", call. = FALSE)
  }

  type <- match.arg(type)

  # error if newdata missing from arguments and object
  if (missing(newdata)) {
    if (is.null(object$newdata)) {
      stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
    } else {
      newdata <- object$newdata
    }
  }

  # deal with newdata_size
  if (missing(newdata_size)) newdata_size <- NULL
  # set newdata_size if needed
  if (is.null(newdata_size) && object$family == "binomial") {
    newdata_size <- rep(1, NROW(newdata))
  }

  # spglm() models a latent Gaussian process w on the link scale via a
  # Laplace approximation (analogous to a GLMM's linear predictor); w plays
  # the role that the observed y plays in conditional.splm() above. Unlike y,
  # w is not observed directly and its own estimation uncertainty is
  # propagated analytically via var_adj below rather than by simulating a new
  # draw of w (see get_conditional_new_from_base_adjust_glm())
  w <- fitted(object, type = "link")
  y <- object$y
  size <- object$size
  base_val_w <- matrix(rep(w, times = samples), ncol = samples)
  if (length(output) == 1 && output == "object") {
    return(base_val_w)
  }
  # handle offset: everything below that is built from Sigma (SqrtSigInv_w,
  # base_val, i.e., the conditional draws) uses the offset-free w, while get_D() below
  # is a derivative of the data model and so must stay on the offset-inclusive
  # linear predictor; see w_offset_free()
  offset_obdata <- model.offset(model.frame(object))
  w <- w_offset_free(w, offset_obdata)


  local_list <- get_local_list_conditional(local, object, newdata)
  if (local_list$approximation == "vecchia" && object$n > 10000) {
    message("local$approximation = \"vecchia\" for spglm() model objects requires a one-time factorization of a dense ", object$n, " x ", object$n, " matrix to account for uncertainty in the latent spatial process. This step does not benefit from vecchia's neighbor truncation (unlike the rest of the simulation) and may be slow and memory-intensive for large observed sample sizes. See Details.")
  }

  betahat <- coef(object)
  X <- model.matrix(object)

  # same composition sampling idea as conditional.splm(): draw beta from its
  # asymptotic sampling distribution N(betahat, vcov(object)) so fixed effect
  # uncertainty propagates into the conditional draws of w below
  cov_betahat_lowchol <- t(chol(vcov(object)))
  new_betahat <- vapply(seq_len(samples), function(x) as.numeric(cov_betahat_lowchol %*% rnorm(length(betahat))), numeric(length(betahat)))
  # beta0 force to matrix
  if (!is.matrix(new_betahat)) {
    new_betahat <- matrix(new_betahat, nrow = 1)
  }
  new_betahat <- sweep(new_betahat, 1, betahat, "+")
  # new_fitted <- X %*% new_betahat
  # new_resid <- w - new_fitted

  # now simulate beta and add
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

  # big data approximation, part 1 (see conditional.splm() for part 2, the
  # newdata blocking, applied identically below): restrict to a
  # spatially-representative base sample, keeping X/w/y/size in sync. Skipped
  # entirely when local_list$approximation == "vecchia", which always conditions on
  # ALL observed data (not subsampled, see get_conditional_vecchia_glm()).
  if (local_list$approximation != "vecchia" && local_list$method_base != "all") {
    object$obdata <- object$obdata[local_list$index$base, , drop = FALSE]
    X <- X[local_list$index$base, , drop = FALSE]
    w <- w[local_list$index$base]
    y <- y[local_list$index$base]
    if (!is.null(size)) {
      size <- size[local_list$index$base]
    }
    # keep the offset aligned with w and y so the get_D() call below can put it back
    if (!is.null(offset_obdata)) {
      offset_obdata <- offset_obdata[local_list$index$base]
    }
  }

  # Covariance components needed by var_adj (applied later, in
  # get_conditional_new_from_base_adjust_glm()/get_conditional_vecchia_glm())
  # -- the analytic adjustment for w's own Laplace-approximate estimation
  # uncertainty:
  #  - SigInv: precision of the spatial covariance matrix of w
  #  - Ptheta: SigInv adjusted for fixed effect estimation uncertainty (the
  #    usual "residual maker" projection SigInv - SigInv X (X'SigInv X)^-1 X'SigInv)
  #  - D: GLM working-weight curvature of the response log-likelihood in w
  #    (see get_D()), i.e. the data's contribution to the Hessian
  #  - cov_lowchol_mH: Cholesky factor of -(D - Ptheta), the negative Hessian
  #    of the joint log-likelihood for w -- var_adj uses its inverse (the
  #    Laplace-approximate posterior covariance of w) to inflate the
  #    predictive variance analytically instead of by simulating a new w
  # This factorization is O(n^3) in whichever data it is computed over -- the
  # (possibly subsampled) base sample for "low-rank", or all observed data for
  # "vecchia", since var_adj reflects
  # uncertainty in the single joint Laplace posterior for w, which is not a
  # spatially-local quantity that vecchia's neighbor truncation can shrink without further investigation
  # (see get_conditional_vecchia_glm() and conditional()'s Details).
  cov_lowchol_base <- t(chol(covmatrix(object)))
  SigInv <- chol2inv(t(cov_lowchol_base))
  SqrtSigInv_X <- forwardsolve(cov_lowchol_base, X)
  SigInv_X <- backsolve(t(cov_lowchol_base), SqrtSigInv_X)
  # solved once and reused for every block's conditional mean below (see
  # get_conditional_new_from_base_adjust_glm()). SqrtSigInv_X
  # above can be used for var_adj
  SqrtSigInv_w <- forwardsolve(cov_lowchol_base, cbind(w))
  cov_betahat <- vcov(object, var_correct = FALSE)
  Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
  # get_D() differentiates the data model, so unlike every spatial quantity
  # above it is evaluated at the offset-inclusive linear predictor
  D <- get_D(
    object$family, if (is.null(offset_obdata)) w else w + as.vector(offset_obdata),
    y, size, as.vector(object$coefficients$dispersion)
  )
  cov_lowchol_mH <- t(chol(Matrix::forceSymmetric(-1 * (D - Ptheta)))) # this is actually the inverse of covariance matrix of w
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)

  # w is held fixed at its fitted (Laplace-mode) value rather than simulated:
  # var_adj (applied in get_conditional_new_from_base_adjust_glm()) already
  # supplies w's own estimation uncertainty analytically, so simulating a new
  # draw of w here on top of that would double-count it. residualizing the
  # (fixed) w against the simulated fixed effect trend plays the same role
  # as new_resid in conditional.splm()
  base_val <- w - X %*% new_betahat

  if (local_list$approximation == "vecchia") {
    # vecchia: every newdata location is simulated sequentially, conditional
    # on ALL observed data plus every earlier-simulated newdata location, with
    # var_adj folded in only between pairs of predicted (never observed)
    # locations (see get_conditional_vecchia_glm())
    new_val <- get_conditional_vecchia_glm(object, newdata, newdata_model, base_val, local_list, samples, SigInv, SigInv_X, wts_beta, cov_lowchol_mH)
  } else {
    if (local_list$method_new != "all") {
      x0 <- lapply(local_list$index$new, function(x) newdata_model[x, , drop = FALSE])
      newdata <- lapply(local_list$index$new, function(x) newdata[x, , drop = FALSE])
    } else {
      x0 <- list(newdata_model)
      newdata <- list(newdata)
    }
    newdata_list <- mapply(x = x0, y = newdata, FUN = function(x, y) list(x0 = x, newdata = y), SIMPLIFY = FALSE)

    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      new_val <- parLapply(cl, newdata_list, get_conditional_new_from_base_adjust_glm, object, SqrtSigInv_w, SqrtSigInv_X, new_betahat, cov_lowchol_base, samples, SigInv, SigInv_X, wts_beta, cov_lowchol_mH)
      cl <- parallel::stopCluster(cl)
    } else {
      new_val <- lapply(newdata_list, get_conditional_new_from_base_adjust_glm, object, SqrtSigInv_w, SqrtSigInv_X, new_betahat, cov_lowchol_base, samples, SigInv, SigInv_X, wts_beta, cov_lowchol_mH)
    }


    new_val <- do.call("rbind", new_val)
    if (local_list$method_new != "all") {
      index_new <- do.call("c", local_list$index$new)
      new_val <- new_val[order(index_new), , drop = FALSE]
    }
  }

  new_val <- newdata_model %*% new_betahat + new_val

  if (!is.null(offset)) {
    new_val <- sweep(new_val, 1, offset, "+")
  }

  # new_val holds link-scale draws up to this point; back-transform to the
  # response scale, or simulate a genuinely new observation, only if asked
  if (type != "link") {
    new_val <- invlink_conditional(new_val, type, dispersion = as.vector(coef(object, type = "dispersion")), family = object$family, newdata_size)
  }

  val <- list(newdata = new_val, beta = new_betahat, object = base_val_w)
  if (length(output) == 1) {
    return(val[[output]])
  } else {
    return(val[output])
  }
}

#' Back-transform conditional simulations from the link scale
#'
#' Converts the link-scale draws returned by the shared code path in
#' \code{conditional.spglm()} to the scale requested via its \code{type}
#' argument: the response-scale mean (\code{type = "response"}), or a newly
#' simulated observation drawn from the response distribution with that mean
#' and the fitted dispersion parameter (\code{type = "new"}).
#'
#' @param mu_link A matrix of link-scale conditional simulations (rows are
#'   \code{newdata} observations, columns are simulation draws).
#' @param type \code{"response"} or \code{"new"} (\code{"link"} is handled by
#'   the caller and never reaches this function).
#' @param dispersion The fitted dispersion parameter from \code{object}.
#' @param family The \code{object} family, e.g. \code{"poisson"}, \code{"Gamma"}.
#' @param newdata_size The binomial size for each \code{newdata} row; only
#'   used when \code{family = "binomial"}.
#'
#' @return A matrix the same shape as \code{mu_link} on the requested scale.
#'
#' @noRd
invlink_conditional <- function(mu_link, type, dispersion, family, newdata_size) {

  # compare to delta method?
  n_newdata <- NROW(mu_link)
  n_sim <- NCOL(mu_link)
  len_sim <- seq(1, n_sim)
  # size = 1 here because binomial rescaling (by the true newdata_size) is
  # applied afterward, once, for both type = "response" and type = "new"
  mu <- invlink(mu_link, family, size = 1)
  rm("mu_link")

  if (type == "response") {
    val <- mu
  } else if (type == "new") {

    # each column of mu is one simulation draw's response-scale mean for
    # every newdata row; for type = "new" a genuinely new observation is
    # drawn per column from the family's response distribution at that mean
    # (and the fitted dispersion), rather than just returning the mean itself
    if (family == "poisson") {
      mu_list <- split(t(mu), len_sim)
      val <- vapply(mu_list, function(x) rpois(n_newdata, x), numeric(n_newdata))
    }

    if (family == "nbinomial") {
      mu_list <- split(t(mu), len_sim)
      val <- vapply(mu_list, function(x) rnbinom(n_newdata, mu = x, size = dispersion), numeric(n_newdata))
    }

    if (family == "Gamma") {
      mu_list <- split(t(mu), len_sim)
      val <- vapply(mu_list, function(x) rgamma(n_newdata, shape = dispersion, scale = x / dispersion), numeric(n_newdata))
    }

    if (family == "inverse.gaussian") {
      mu_list <- split(t(mu), len_sim)
      val <- vapply(mu_list, function(x) {
        dispersion_true <- 1 / (x * dispersion)
        statmod::rinvgauss(n_newdata, mean = x, dispersion = dispersion_true)
      }, numeric(n_newdata))
    }

    if (family == "binomial") {
      mu_list <- split(t(mu), len_sim)
      val <- vapply(mu_list, function(x) rbinom(n_newdata, newdata_size, x), numeric(n_newdata))
    }

    if (family == "beta") {
      mu_list <- split(t(mu), len_sim)
      val <- vapply(mu_list, function(x) {
        a <- x * dispersion
        b <- (1 - x) * dispersion
        val <- rbeta(n_newdata, shape1 = a, shape2 = b)
        val <- pmax(1e-4, val)
        val <- pmin(1 - 1e-4, val)
      }, numeric(n_newdata))
    }
  }

  # invlink() above returned a proportion (size = 1); rescale to counts using
  # the actual newdata_size now that both type branches have produced their
  # response-scale value (rbinom() already returns counts, so "new" skips this)
  if (type == "response" && family == "binomial") {
    val <- sweep(val, 1, newdata_size, "*")
  }

  val
}
