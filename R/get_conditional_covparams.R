#' Retrieve the covariance matrix of the covariance parameter estimates for
#' \code{conditional()}'s \code{simulate_covparams = TRUE} path
#'
#' Requires \code{object$vcov$cov} (the same covariance matrix
#' \code{satterthwaite()} uses, populated at fit time whenever \code{ddf =
#' "satterthwaite"} -- see \code{get_fit_ddf()}/\code{satterthwaite_core()}) to
#' already exist, rather than ever recomputing it here.
#'
#' @param object A fitted \code{splm} model object.
#'
#' @return \code{object$vcov$cov}, or an informative error if it is \code{NULL}.
#'
#' @noRd
get_vcov_theta_for_conditional <- function(object) {
  vcov_theta <- object$vcov$cov
  if (is.null(vcov_theta)) {
    # NOTE: could recompute this on the fly via satterthwaite_core() (the way
    # get_emmeans_dffun() does), but that repeats an expensive Hessian-based
    # computation on top of what simulate_covparams's own per-draw cost already
    # warns about.
    stop("simulate_covparams = TRUE requires object$vcov$cov, the covariance matrix of the covariance parameter estimates. This is only available when the model was fit with ddf = \"satterthwaite\" (the default for n <= 500; see splm()'s ddf argument). Refit with ddf = \"satterthwaite\" to use simulate_covparams = TRUE.", call. = FALSE)
  }
  vcov_theta
}

#' Substitute a simulated free-parameter draw into a covariance parameter
#' object, validated
#'
#' Reuses \code{spcov_params()}'s own \code{stop()}-based validation (every
#' family/bound check it already performs) as the validity gate for a
#' simulated draw, instead of duplicating that constraint logic.
#' \code{randcov_params()} does not itself validate non-negativity, so that
#' check is added explicitly.
#'
#' @param draw A named numeric vector (names a subset of
#'   \code{c(spcov_names_free, randcov_names_free)}) of simulated values for
#'   the free covariance parameters.
#' @param spcov_type The fitted \code{spcov_type} (e.g. \code{"exponential"}).
#' @param spcov_names_free,randcov_names_free Character vectors of free
#'   parameter names (random effect names empty/character(0) if none).
#' @param spcov_params_full,randcov_params_full The fitted parameter vectors
#'   (free and known), providing values for any known/fixed entries.
#'
#' @return \code{list(spcov_params_b, randcov_params_b)} if \code{draw} is
#'   valid, or \code{NULL} if any constraint is violated (signals "redraw").
#'
#' @noRd
try_build_theta <- function(draw, spcov_type, spcov_names_free, randcov_names_free,
                             spcov_params_full, randcov_params_full) {
  spcov_params_new <- spcov_params_full
  if (length(spcov_names_free)) {
    spcov_params_new[spcov_names_free] <- draw[spcov_names_free]
  }
  randcov_params_new <- randcov_params_full
  if (length(randcov_names_free)) {
    randcov_params_new[randcov_names_free] <- draw[randcov_names_free]
  }

  spcov_params_b <- tryCatch(
    do.call(spcov_params, c(list(spcov_type = spcov_type), as.list(spcov_params_new))),
    error = function(e) NULL
  )
  if (is.null(spcov_params_b)) {
    return(NULL)
  }

  if (!is.null(randcov_params_new)) {
    if (any(randcov_params_new < 0)) {
      return(NULL)
    }
    randcov_params_b <- do.call(randcov_params, as.list(randcov_params_new))
  } else {
    randcov_params_b <- NULL
  }

  list(spcov_params_b = spcov_params_b, randcov_params_b = randcov_params_b)
}

#' Clamp a simulated free-parameter draw to the nearest valid boundary
#'
#' Fallback used once \code{simulate_theta_draw()}'s reject-and-redraw attempts
#' are exhausted. This mirrors \code{spcov_params()}'s own bounds so the result is
#' guaranteed valid by construction. A clamped \code{range} of exactly 0 also
#' forces \code{de} to 0, since a zero-range spatial term is a degenerate
#' pure-nugget limit and leaving \code{de > 0} there risks division by zero in
#' \code{exp(-d/range)}-style formulas.
#'
#' @inheritParams try_build_theta
#'
#' @return \code{list(spcov_params_b, randcov_params_b)}, always valid.
#'
#' @noRd
clamp_theta_draw <- function(draw, spcov_type, spcov_names_free, randcov_names_free,
                              spcov_params_full, randcov_params_full) {
  spcov_params_new <- spcov_params_full
  if (length(spcov_names_free)) {
    spcov_params_new[spcov_names_free] <- draw[spcov_names_free]
  }
  randcov_params_new <- randcov_params_full
  if (length(randcov_names_free)) {
    randcov_params_new[randcov_names_free] <- draw[randcov_names_free]
  }

  if ("de" %in% names(spcov_params_new)) {
    spcov_params_new["de"] <- max(0, spcov_params_new["de"])
  }
  if ("ie" %in% names(spcov_params_new)) {
    spcov_params_new["ie"] <- max(0, spcov_params_new["ie"])
  }
  if ("range" %in% names(spcov_params_new) && !spcov_type %in% c("none", "ie")) {
    spcov_params_new["range"] <- max(0, spcov_params_new["range"])
    if (spcov_params_new["range"] == 0) {
      spcov_params_new["de"] <- 0
    }
  }
  if ("rotate" %in% names(spcov_params_new)) {
    spcov_params_new["rotate"] <- min(max(spcov_params_new["rotate"], 0), pi)
  }
  if ("scale" %in% names(spcov_params_new)) {
    spcov_params_new["scale"] <- min(max(spcov_params_new["scale"], 0), 1)
  }
  if ("extra" %in% names(spcov_params_new)) {
    extra_val <- spcov_params_new["extra"]
    if (spcov_type == "matern") {
      extra_val <- min(max(extra_val, 0.2), 5)
    } else if (spcov_type == "cauchy") {
      extra_val <- max(extra_val, .Machine$double.eps)
    } else if (spcov_type == "pexponential") {
      extra_val <- min(max(extra_val, .Machine$double.eps), 2)
    }
    spcov_params_new["extra"] <- extra_val
  }

  spcov_params_b <- do.call(spcov_params, c(list(spcov_type = spcov_type), as.list(spcov_params_new)))

  if (!is.null(randcov_params_new)) {
    randcov_params_new <- pmax(0, randcov_params_new)
    randcov_params_b <- do.call(randcov_params, as.list(randcov_params_new))
  } else {
    randcov_params_b <- NULL
  }

  list(spcov_params_b = spcov_params_b, randcov_params_b = randcov_params_b)
}

#' Simulate one valid draw of the free covariance parameters
#'
#' Reject-and-redraw up to \code{max_attempts} times (via
#' \code{try_build_theta()}'s reuse of \code{spcov_params()}'s own
#' validation); if every attempt is invalid, falls back to
#' \code{clamp_theta_draw()} so this always returns a usable value.
#'
#' @inheritParams try_build_theta
#' @param theta_hat_free A named numeric vector: the fitted values of the free
#'   covariance parameters (names \code{c(spcov_names_free,
#'   randcov_names_free)}).
#' @param vcov_theta_lowchol The lower-triangular Cholesky factor of
#'   \code{theta_hat_free}'s covariance matrix.
#' @param max_attempts The number of reject-and-redraw attempts before falling
#'   back to clamping. The default is 50.
#'
#' @return \code{list(spcov_params_b, randcov_params_b)}.
#'
#' @noRd
simulate_theta_draw <- function(theta_hat_free, vcov_theta_lowchol, spcov_type,
                                 spcov_names_free, randcov_names_free,
                                 spcov_params_full, randcov_params_full, max_attempts = 50) {
  draw <- NULL
  for (attempt in seq_len(max_attempts)) {
    draw <- theta_hat_free + as.numeric(vcov_theta_lowchol %*% rnorm(length(theta_hat_free)))
    out <- try_build_theta(
      draw, spcov_type, spcov_names_free, randcov_names_free,
      spcov_params_full, randcov_params_full
    )
    if (!is.null(out)) {
      return(out)
    }
  }
  clamp_theta_draw(
    draw, spcov_type, spcov_names_free, randcov_names_free,
    spcov_params_full, randcov_params_full
  )
}

#' Conditionally simulate \code{newdata}/beta draws while also propagating
#' covariance parameter estimation uncertainty
#'
#' Implements \code{conditional.splm()}'s \code{simulate_covparams = TRUE}
#' path: for each of \code{samples} draws, simulate a new covariance parameter
#' vector (\code{\link{simulate_theta_draw}}), recompute the GLS covariance of
#' betahat under it, draw a new betahat, and redo the conditional spatial draw
#' with that draw's own covariance parameters. Unlike the rest of
#' \code{conditional()}, this cannot vectorize across samples (nothing can be
#' shared once the covariance parameters change every draw), so it is a plain
#' sequential loop, one fresh \code{O(n^3)} factorization per sample.
#'
#' Only called when \code{local_list$method_base == local_list$method_new ==
#' "all"} is guaranteed (see \code{conditional.splm()}'s \code{local_active}
#' check), so this always operates on the full observed data and all of
#' \code{newdata} at once and there is no base-subsampling/blocking to account for.
#'
#' @param object A fitted \code{splm} model object.
#' @param newdata The (already processed, i.e. post-
#'   \code{get_newdata_model_matrix()}) \code{newdata} data frame.
#' @param y The (offset-adjusted) observed response.
#' @param X The observed data design matrix.
#' @param betahat_hat The fitted coefficient vector (the center every drawn
#'   betahat is simulated around, regardless of the draw's own covariance
#'   parameters).
#' @param vcov_theta The covariance matrix of the free covariance parameter
#'   estimates (\code{object$vcov$cov}, see
#'   \code{\link{get_vcov_theta_for_conditional}}).
#' @param samples The number of simulations.
#'
#' @return A list with elements \code{new_betahat} (a p x \code{samples}
#'   matrix), \code{new_val} (a \code{NROW(newdata)} x \code{samples} matrix
#'   of mean-zero conditional spatial residual draws, matching the shape
#'   \code{conditional.splm()}'s existing paths return before the fixed
#'   effect trend is added back in), and the simulated covariance parameter
#'   draws themselves -- \code{new_cov} (all of them, spcov then randcov
#'   rows, mirroring \code{vcov(object, type = "cov")}), \code{new_spcov},
#'   and \code{new_randcov} (\code{NULL} if the model has no random effects)
#'   -- each a (parameter) x \code{samples} matrix.
#'
#' @noRd
get_conditional_covparams <- function(object, newdata, y, X, betahat_hat, vcov_theta, samples) {

  spcov_params_full <- coef(object, type = "spcov")
  spcov_type <- class(spcov_params_full)
  spcov_is_known <- object$is_known$spcov
  spcov_names_free <- names(spcov_is_known)[!spcov_is_known]

  has_randcov <- !is.null(object$random)
  if (has_randcov) {
    randcov_params_full <- coef(object, type = "randcov")
    randcov_is_known <- object$is_known$randcov
    randcov_names_free <- names(randcov_is_known)[!randcov_is_known]
  } else {
    randcov_params_full <- NULL
    randcov_names_free <- character(0)
  }

  cov_names_free <- c(spcov_names_free, randcov_names_free)
  cov_val_full <- c(as.numeric(spcov_params_full), as.numeric(randcov_params_full))
  names(cov_val_full) <- c(names(spcov_params_full), names(randcov_params_full))
  theta_hat_free <- cov_val_full[cov_names_free]

  # order/subset to match cov_names_free exactly, mirroring
  # get_cov_gradients_context.splm()'s own derivation of these same names
  vcov_theta <- vcov_theta[cov_names_free, cov_names_free, drop = FALSE]
  vcov_theta_lowchol <- t(chol(vcov_theta))

  newdata_n <- NROW(newdata)
  p <- length(betahat_hat)
  new_betahat <- matrix(NA_real_, p, samples)
  new_val <- matrix(NA_real_, newdata_n, samples)
  new_spcov <- matrix(NA_real_, length(spcov_params_full), samples, dimnames = list(names(spcov_params_full), NULL))
  new_randcov <- if (has_randcov) {
    matrix(NA_real_, length(randcov_params_full), samples, dimnames = list(names(randcov_params_full), NULL))
  } else {
    NULL
  }

  for (b in seq_len(samples)) {
    theta_b <- simulate_theta_draw(
      theta_hat_free, vcov_theta_lowchol, spcov_type,
      spcov_names_free, randcov_names_free, spcov_params_full, randcov_params_full
    )

    object_b <- object
    object_b$coefficients$spcov <- theta_b$spcov_params_b
    object_b$coefficients$randcov <- theta_b$randcov_params_b

    # pure nugget (independent error, no spatial dependence or random
    # effects): the base covariance matrix is diagonal, so a plain sqrt()
    # gives its Cholesky factor without paying for a full chol() -- same
    # optimization as the existing (simulate_covparams = FALSE) code path
    pure_nugget <- theta_b$spcov_params_b[["de"]] == 0 && is.null(theta_b$randcov_params_b)
    if (pure_nugget) {
      cov_lowchol_base_b <- Matrix::Diagonal(n = object$n, x = sqrt(theta_b$spcov_params_b[["ie"]]))
    } else {
      cov_lowchol_base_b <- t(chol(covmatrix(object_b)))
    }

    # recompute betahat's GLS covariance under theta_b and draw a new betahat.
    # REML and ML share this same formula given theta, so no estmethod
    # branch is needed here
    SqrtSigInv_X_b <- forwardsolve(cov_lowchol_base_b, X)
    cov_betahat_b <- chol2inv(chol(crossprod(SqrtSigInv_X_b)))
    betahat_b <- betahat_hat + as.numeric(t(chol(cov_betahat_b)) %*% rnorm(p))
    resid_b <- y - X %*% betahat_b

    # conditional draw under theta_b for this one sample
    cov_base_new_b <- covmatrix(object_b, newdata, cov_type = "obs.pred")
    cov_new_b <- covmatrix(object_b, newdata, cov_type = "pred.pred")
    SqrtSigInv_c0_b <- forwardsolve(cov_lowchol_base_b, cov_base_new_b)
    cond_cov_b <- cov_new_b - crossprod(SqrtSigInv_c0_b)
    chol_cond_cov_b <- t(chol(cond_cov_b))
    cond_mu_b <- crossprod(SqrtSigInv_c0_b, forwardsolve(cov_lowchol_base_b, resid_b))

    new_betahat[, b] <- betahat_b
    new_val[, b] <- as.numeric(chol_cond_cov_b %*% rnorm(newdata_n)) + as.numeric(cond_mu_b)
    new_spcov[, b] <- as.numeric(theta_b$spcov_params_b)
    if (has_randcov) {
      new_randcov[, b] <- as.numeric(theta_b$randcov_params_b)
    }
  }

  # "cov" mirrors vcov(object, type = "cov")'s own joint spcov+randcov
  # convention -- the same row order get_cov_gradients_context.splm() uses
  # (spcov names, then randcov names)
  new_cov <- rbind(new_spcov, new_randcov)

  list(new_betahat = new_betahat, new_val = new_val, new_cov = new_cov, new_spcov = new_spcov, new_randcov = new_randcov)
}
