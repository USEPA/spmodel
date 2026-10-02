#' @rdname conditional
#' @method conditional spautor
#' @export
conditional.spautor <- function(object, newdata, output = "newdata", samples = 1000, ...) {
  output <- check_conditional_areal(object, output, samples, list(...))
  context <- get_prediction_object_spautor(object, newdata, scale = NULL, local = FALSE)
  y <- model.response(model.frame(object))
  if (identical(output, "object")) {
    return(conditional_areal_snapshot(object, y, samples))
  }

  cov_context <- get_conditional_areal_cov(object)
  betahat <- coef(object)
  new_betahat <- t(chol(vcov(object))) %*%
    matrix(rnorm(length(betahat) * samples), length(betahat), samples)
  new_betahat <- sweep(new_betahat, 1, betahat, "+")
  rownames(new_betahat) <- names(betahat)
  residual <- w_offset_free(y, model.offset(model.frame(object))) -
    model.matrix(object) %*% new_betahat
  new_val <- draw_conditional_areal(context, cov_context, residual, new_betahat, samples)
  conditional_areal_output(object, new_val, new_betahat, y, output, samples)
}

#' Validate exact areal simulation settings
#' @noRd
check_conditional_areal <- function(object, output, samples, dots) {
  if ("local" %in% names(dots)) {
    stop("local is not supported for autoregressive conditional simulation; only exact simulation is available.", call. = FALSE)
  }
  if ("simulate_covparams" %in% names(dots)) {
    stop("simulate_covparams is not supported for autoregressive conditional simulation.", call. = FALSE)
  }
  if (length(dots)) stop("Unused arguments in autoregressive conditional simulation.", call. = FALSE)
  if (!is.numeric(samples) || length(samples) != 1L || !is.finite(samples) ||
      samples < 1 || samples != floor(samples)) {
    stop("samples must be a positive integer.", call. = FALSE)
  }
  if (!is.character(output) || !length(output) || anyNA(output) ||
      any(!output %in% c("newdata", "beta", "object", "all"))) {
    stop('output must be "newdata", "beta", "object", or "all".', call. = FALSE)
  }
  if ("all" %in% output) output <- c("newdata", "beta", "object")
  if (is.null(object$newdata) || !NROW(object$newdata) || !length(object$missing_index)) {
    stop("No missing data to simulate. Fit the model with NA response values for the locations you want to simulate.", call. = FALSE)
  }
  output
}

#' Prepare observed and missing covariance factors without changing the graph
#' @noRd
get_conditional_areal_cov <- function(object) {
  randcov_Zs <- get_randcov_Zs(object$data, get_randcov_names(object$random))
  cov_full <- as.matrix(cov_matrix(coef(object, type = "spcov"), object$W,
    coef(object, type = "randcov"), randcov_Zs,
    partition_matrix(object$partition_factor, object$data), object$M))
  observed <- object$observed_index
  missing <- object$missing_index
  cov_lowchol <- t(chol(cov_full[observed, observed, drop = FALSE]))
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol, cov_full[observed, missing, drop = FALSE])
  cond_cov <- cov_full[missing, missing, drop = FALSE] - crossprod(SqrtSigInv_c0)
  list(cov_lowchol = cov_lowchol, SqrtSigInv_c0 = SqrtSigInv_c0,
    cond_lowchol = t(chol(cond_cov)))
}

#' Draw all missing-site residuals jointly across simulation columns
#' @noRd
draw_conditional_areal <- function(context, cov_context, residual, new_betahat, samples,
                                   wts_residual = NULL) {
  n_new <- NROW(context$newdata_model)
  spatial_mean <- if (is.null(wts_residual)) {
    crossprod(cov_context$SqrtSigInv_c0, forwardsolve(cov_context$cov_lowchol, residual))
  } else crossprod(wts_residual, residual)
  new_val <- context$newdata_model %*% new_betahat +
    spatial_mean +
    cov_context$cond_lowchol %*% matrix(rnorm(n_new * samples), n_new, samples)
  if (!is.null(context$offset)) new_val <- sweep(new_val, 1, context$offset, "+")
  new_val
}

#' Replicate observed values in fitted row order
#' @noRd
conditional_areal_snapshot <- function(object, observed, samples) {
  matrix(rep(observed, samples), nrow = length(observed), ncol = samples,
    dimnames = list(as.character(object$observed_index), NULL))
}

#' Assemble areal simulation outputs
#' @noRd
conditional_areal_output <- function(object, new_val, new_betahat, observed, output, samples) {
  new_val <- matrix(new_val, nrow = length(object$missing_index), ncol = samples,
    dimnames = list(as.character(object$missing_index), NULL))
  val <- list(newdata = new_val, beta = new_betahat)
  if ("object" %in% output) val$object <- conditional_areal_snapshot(object, observed, samples)
  if (length(output) == 1L) val[[output]] else val[output]
}
