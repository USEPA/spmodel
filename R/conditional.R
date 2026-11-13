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
#'   \code{c("newdata", "beta", "object")} (\code{"all"} is shorthand for this full vector).
#'   If \code{simulate_covparams = TRUE},
#'   the output type can also be any subset of \code{c("cov", "spcov", "randcov")}.
#'   For generalized linear models, \code{"latent"} returns simulated observed-site
#'   link values, including offsets.
#'   \code{"random"} returns realized random effects.
#'    The default is \code{"newdata"}. See Details for more.
#' @param type For \code{spglm()} and \code{spgautor()} model objects, the scale of the conditional
#'   simulations for \code{newdata}.
#'   When \code{type = "link"}, the predicted means on the
#'   link scale are returned. When \code{type = "response"}, the predicted means
#'   on the response scale are returned. When \code{type = "new"}, a new observation
#'   is simulated from the appropriate response distribution with mean equal to
#'   the mean on the response scale and dispersion equal to the dispersion parameter
#'   from \code{object}. The default is \code{"link"}.
#' @param samples  The number of conditional simulations. The default is
#'   \code{1,000}.
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
#'   observed responses (\code{splm()} or \code{spautor()}) or fitted link values (including offsets)
#'   (\code{spglm()} or \code{spgautor()}) are returned in observed-row order,
#'   repeated across simulation columns. If \code{"latent"} is requested, the
#'   actual observed link-scale draws are returned, with offsets included.
#'  \code{type} affects only \code{newdata}.
#'   For example, \code{c("newdata", "beta")} returns
#'   the conditional simulations both for \code{newdata} and for the
#'   fixed effects. If \code{"cov"}/\code{"spcov"}/\code{"randcov"} is in
#'   \code{output} (only available when \code{simulate_covparams = TRUE}),
#'   the simulated covariance parameter draws themselves are returned.
#'
#' @return If \code{output = "newdata"}, an a x b matrix of conditional simulations
#'   for each row in \code{newdata}, where a is the
#'   number of rows in \code{newdata} and b is the number of samples.
#'   If \code{output = "beta"}, an p x b matrix of conditional simulations for each
#'   element in \code{coef(object)}, where p is the
#'   number of fixed effects and b is the number of samples.
#'   If \code{output = "object"}, an n x b matrix of observed responses
#'   (\code{splm()} and \code{spautor()}) or fitted link values (including offsets)
#'   (\code{spglm()} and \code{spgautor()}), where n is the number of observed rows and b is the number of samples.
#'   \code{output = "latent"} returns an n x b matrix of simulated observed
#'   link values (including offsets).
#'   \code{output = "random"} returns a q x b matrix of realized random effects coefficients,
#'   separated by partition factor (when relevant).
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
conditional.splm <- function(object, newdata, output = "newdata", samples = 1000, simulate_covparams = FALSE, ...) {

  if (!is.numeric(samples) || length(samples) != 1L || !is.finite(samples) || samples < 1 || samples != floor(samples)) {
    stop("samples must be a positive integer.", call. = FALSE)
  }

  if (!is.logical(simulate_covparams) || length(simulate_covparams) != 1 || is.na(simulate_covparams)) {
    stop("simulate_covparams must be TRUE or FALSE.", call. = FALSE)
  }

  if ("all" %in% output) {
    output <- c("newdata", "beta", "object", if ("random" %in% output) "random")
  }
  if (!is.character(output) || !length(output) || anyNA(output) || any(!output %in% c("newdata", "beta", "object", "cov", "spcov", "randcov", "random"))) {
    stop("output must be \"newdata\", \"beta\", \"object\", \"cov\", \"spcov\", \"randcov\", \"random\", or \"all\".", call. = FALSE)
  }
  if ("random" %in% output && is.null(object$random)) {
    stop('output = "random" requires a model with random effects.', call. = FALSE)
  }

  # error if newdata missing from arguments and object
  if (missing(newdata)) {
    if (!"newdata" %in% output) {
      newdata <- object$obdata[FALSE, , drop = FALSE]
    } else if (is.null(object$newdata)) {
      stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
    } else {
      newdata <- object$newdata
    }
  }

  coordinates <- prepare_conditional_coordinates(object, newdata)
  object <- coordinates$object
  newdata <- coordinates$newdata
  iid_variance <- get_conditional_iid_variance(object, simulate_covparams)
  design <- if (!is.null(object$random)) get_conditional_random_design(object, newdata) else NULL
  random_val <- NULL

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
  base_val_y <- if ("object" %in% output) matrix(rep(y, times = samples), ncol = samples) else NULL
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
  newdata_model_list <- if (NROW(newdata)) get_newdata_model_matrix(object, newdata) else
    list(newdata = object$obdata[FALSE, , drop = FALSE], newdata_model = X[FALSE, , drop = FALSE], offset = NULL)
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
  if (!is.null(iid_variance) && !is.null(object$partition_factor) && NROW(newdata)) {
    model.frame(object$partition_factor, newdata)
  }

  if (isTRUE(simulate_covparams)) {
    # covariance parameters resimulated per draw -- see get_conditional_covparams()
    # for why this cannot share the vectorized-across-samples approach below
    covparams_out <- get_conditional_covparams(object, newdata, y, X, betahat, vcov_theta, samples, design, "random" %in% output)
    new_betahat <- covparams_out$new_betahat
    new_val <- covparams_out$new_val
    new_cov <- covparams_out$new_cov
    new_spcov <- covparams_out$new_spcov
    new_randcov <- covparams_out$new_randcov
    if ("random" %in% output) random_val <- conditional_random_output(covparams_out$random, design)
  } else if (!is.null(object$random)) {
    draws <- get_conditional_random_lm(object, newdata, newdata_model, y, samples, design, "random" %in% output)
    new_betahat <- draws$beta
    new_val <- draws$newdata
    random_val <- draws$random
    new_cov <- new_spcov <- new_randcov <- NULL
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
    chunk_size <- samples

    if (!NROW(newdata)) {
      new_val <- matrix(numeric(), 0, samples)
    } else if (!is.null(iid_variance)) {
      new_val <- draw_conditional_iid_error(NROW(newdata), samples, iid_variance)
    } else {
      cov_lowchol_base <- t(chol(covmatrix(object)))
      SqrtSigInv_y <- forwardsolve(cov_lowchol_base, cbind(y))
      SqrtSigInv_X <- forwardsolve(cov_lowchol_base, X)
      new_val <- draw_conditional_exact_lm(
        newdata, object, SqrtSigInv_y, SqrtSigInv_X,
        new_betahat, cov_lowchol_base, samples, chunk_size)
    }
  }

  # add the simulated fixed effect trend (X0 %*% beta_b) back onto the
  # (mean-zero) conditional spatial residual draws -- completes the
  # composition sampling described previously
  if (isTRUE(simulate_covparams) || is.null(object$random)) {
    chunk_size <- if (!is.null(iid_variance)) 1000L else samples
    for (start in seq.int(1L, samples, by = chunk_size)) {
      cols <- seq.int(start, min(samples, start + chunk_size - 1L))
      new_val[, cols] <- newdata_model %*% new_betahat[, cols, drop = FALSE] + new_val[, cols, drop = FALSE]
    }
  }

  if (!is.null(offset)) {
    new_val <- sweep(new_val, 1, offset, "+")
  }

  if ("random" %in% output) random_val <- check_conditional_random_output(random_val)
  val <- list(newdata = new_val, beta = new_betahat, object = base_val_y, cov = new_cov, spcov = new_spcov, randcov = new_randcov, random = random_val)
  if (length(output) == 1) {
    return(val[[output]])
  } else {
    return(val[output])
  }
}

#' @rdname conditional
#' @method conditional spglm
#' @export
conditional.spglm <- function(object, newdata, output = "newdata", type = c("link", "response", "new"), samples = 1000, newdata_size, ...) {

  if (!is.numeric(samples) || length(samples) != 1L || !is.finite(samples) || samples < 1 || samples != floor(samples)) {
    stop("samples must be a positive integer.", call. = FALSE)
  }


  if ("all" %in% output) {
    output <- c("newdata", "beta", "object", if ("latent" %in% output) "latent", if ("random" %in% output) "random")
  }
  if (!is.character(output) || !length(output) || anyNA(output) || any(!output %in% c("newdata", "beta", "object", "latent", "random"))) {
    stop("output must be \"newdata\", \"beta\", \"object\", \"latent\", \"random\", or \"all\".", call. = FALSE)
  }
  if ("random" %in% output && is.null(object$random)) {
    stop('output = "random" requires a model with random effects.', call. = FALSE)
  }

  type <- match.arg(type)

  # Latent-only calls need no prediction data; explicit prediction calls retain
  # the existing missing-response/newdata behavior.
  if (missing(newdata)) {
    if (!"newdata" %in% output) {
      newdata <- NULL
    } else if (is.null(object$newdata)) {
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

  base_val_w <- if ("object" %in% output) {
    matrix(rep(fitted(object, type = "link"), times = samples), ncol = samples)
  } else NULL
  if (length(output) == 1 && output == "object") {
    return(base_val_w)
  }
  iid_variance <- get_conditional_iid_variance(object)
  coordinates <- prepare_conditional_coordinates(object, newdata)
  object <- coordinates$object
  newdata <- coordinates$newdata
  design <- if (!is.null(object$random)) get_conditional_random_design(object, newdata) else NULL
  random_val <- NULL
  if (is.null(iid_variance) && object$n > 10000) {
    message("Exact spatial generalized linear model conditional simulation requires a dense observed latent precision factorization. See Details.")
  }
  newdata_model_list <- if (NROW(newdata)) get_newdata_model_matrix(object, newdata) else
    list(newdata = object$obdata[FALSE, , drop = FALSE], newdata_model = model.matrix(object)[FALSE, , drop = FALSE], offset = NULL)
  newdata <- newdata_model_list$newdata
  newdata_model <- newdata_model_list$newdata_model
  offset <- newdata_model_list$offset
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts
  if (!is.null(iid_variance) && !is.null(object$partition_factor) && NROW(newdata)) {
    model.frame(object$partition_factor, newdata)
  }

  if (!is.null(iid_variance)) {
    draws <- get_conditional_iid_glm(object, newdata_model, iid_variance, samples, latent = "latent" %in% output)
    new_val <- draws$newdata
    new_betahat <- draws$beta
  } else {
    joint <- get_conditional_glm_joint(object)
    draws <- draw_conditional_glm_joint(joint, samples, residual = TRUE)
    new_betahat <- draws$beta
    if (!is.null(object$random)) {
      effect <- prepare_conditional_random(object, design)
      prediction <- prepare_conditional_random_prediction(effect, newdata)
      a <- draw_conditional_random(effect, draws$residual)
      new_val <- newdata_model %*% new_betahat + draw_conditional_random_prediction(effect, prediction, draws$residual, a)
      if ("random" %in% output) random_val <- conditional_random_output(a, design)
    } else new_val <- if (NROW(newdata)) newdata_model %*% new_betahat +
      draw_conditional_exact_glm(newdata, object, draws$residual, joint$cov_lowchol, samples) else matrix(numeric(), 0, samples)
    if ("latent" %in% output) draws$latent <- draws$residual + model.matrix(object) %*% new_betahat
    rm(joint)
  }

  # The latent output is on the link scale but needs offsets added back.
  latent_val <- if ("latent" %in% output) draws$latent else NULL
  if (!is.null(latent_val)) {
    observed_offset <- model.offset(model.frame(object))
    if (!is.null(observed_offset)) latent_val <- sweep(latent_val, 1L, observed_offset, "+")
    rownames(latent_val) <- rownames(model.matrix(object))
  }

  if (!is.null(offset)) {
    new_val <- sweep(new_val, 1, offset, "+")
  }

  # new_val holds link-scale draws up to this point; back-transform to the
  # response scale, or simulate a genuinely new observation, only if asked
  if (type != "link" && NROW(new_val)) {
    chunk_size <- samples
    for (start in seq.int(1L, samples, by = chunk_size)) {
      cols <- seq.int(start, min(samples, start + chunk_size - 1L))
      new_val[, cols] <- invlink_conditional(new_val[, cols, drop = FALSE], type,
        dispersion = as.vector(coef(object, type = "dispersion")), family = object$family, newdata_size)
    }
  }

  if ("random" %in% output) random_val <- check_conditional_random_output(random_val)
  val <- list(newdata = new_val, beta = new_betahat, object = base_val_w, latent = latent_val, random = random_val)
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

  if (is.null(dim(val))) dim(val) <- c(n_newdata, n_sim)
  val
}

prepare_conditional_coordinates <- function(object, newdata) {
  if (inherits(newdata, "sf")) {
    newdata <- sf_to_df(suppressWarnings(sf::st_centroid(newdata)))
    names(newdata)[names(newdata) == ".xcoord"] <- object$xcoord
    names(newdata)[names(newdata) == ".ycoord"] <- object$ycoord
  }
  if (!inherits(coef(object, type = "spcov"), c("none", "ie"))) {
    if (object$dim_coords == 1) {
      # Keep ignored coordinates available as predictors while using 1-D distances.
      coordinate_names <- make.unique(c(names(object$obdata), names(newdata), ".conditional_ycoord"))
      object$ycoord <- coordinate_names[[length(coordinate_names)]]
      object$obdata[[object$ycoord]] <- rep(0, NROW(object$obdata))
      newdata[[object$ycoord]] <- rep(0, NROW(newdata))
    }
    if (NROW(newdata)) check_newdata_coords(newdata, object$xcoord, object$ycoord)
  }
  list(object = object, newdata = newdata)
}
