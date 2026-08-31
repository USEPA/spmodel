#' Predict each model term's contribution to the linear predictor
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()]
#' @param X_newdata The design matrix for the data requiring prediction
#' @param se.fit Whether to compute standard errors for each term
#' @param scale A scale multiplier for the standard errors (or \code{NULL})
#' @param df Degrees of freedom used for interval quantiles
#' @param interval The type of interval (\code{"none"} or \code{"confidence"}; \code{"prediction"} is not supported)
#' @param level The confidence interval level
#' @param add_newdata_rows Whether to name output rows using \code{object$missing_index}
#'   rather than \code{rownames(X_newdata)}
#' @param terms Optional subset of terms to return
#' @param ... Additional arguments passed to \code{vcov()}
#'
#' @return A matrix (or list of matrices, if \code{se.fit} or \code{interval == "confidence"})
#'   with one column per model term (collapsing over factor levels within a
#'   term) and one row per prediction, mirroring \code{predict.lm(type = "terms")}
#'
#' @noRd
predict_terms <- function(object, X_newdata, se.fit, scale, df, interval, level, add_newdata_rows, terms, ...) {
  if (interval == "prediction") {
    stop('If type = "terms", interval must be "none" or "confidence".', call. = FALSE)
  }

  # get the design matrix for the data
  X <- model.matrix(object)
  # get the column means for the observed data
  avx <- colMeans(X)
  # get the regression coefficients
  beta <- matrix(coef(object), ncol = 1)
  # get the covariance matrix of the regression coefficients
  vc <- vcov(object, ...)
  # the predict.lm function sweeps the newdata by the means of the observed data
  # that is, it subracts the column means of observed data from the columns of
  # the prediction data
  hasIntercept <- attr(terms(object), "intercept") > 0L
  if (hasIntercept) {
    X_newdata_cent <- sweep(X_newdata, 2L, avx)
    # the intercept needs to be adjusted due to centering of covariates
    constant <- as.numeric(crossprod(avx, beta))
  } else {
    X_newdata_cent <- X_newdata
    constant <- 0
  }
  # the number of terms (not including intercept)
  assign <- attr(X, "assign")
  nterms <- max(assign)
  # if there is an intercept (you will need to create code to test and adjust
  # for the no intercept case
  fit <- matrix(numeric(0), nrow = NROW(X_newdata), ncol = nterms)
  colnames(fit) <- attr(terms(object), "term.labels")
  if (add_newdata_rows) {
    rownames(fit) <- object$missing_index
  } else {
    rownames(fit) <- rownames(X_newdata)
  }

  if (se.fit || interval == "confidence") {
    se <- fit
  }
  # the number of columns is the number of terms in the formula, so it collapses
  # over factor levels
  for (i in seq_len(nterms)) {
    X_index <- attr(X, "assign") == i
    # the fits per row are just the centered x-values for each term
    # times the regression coefficients, and then summed if these are factors
    fit[, i] <- X_newdata_cent[, X_index, drop = FALSE] %*% beta[X_index, , drop = FALSE]
    # get standard errors
    if (se.fit || interval == "confidence") {
      X_newdata_cent_sub <- X_newdata_cent[, X_index, drop = FALSE]
      vc_sub <- vc[X_index, X_index, drop = FALSE]
      # the fits are just linear combinations, so standard variance rules apply
      se[, i] <- sqrt(get_diag_XVXt(X_newdata_cent_sub, vc_sub))
    }
  }

  fit <- structure(fit, constant = constant)
  out <- list(fit = fit)
  if (se.fit || interval == "confidence") {
    # for some reason, predict.lm returns se even if se.fit = FALSE when
    # interval is confidence
    out$se.fit <- se
    if (!is.null(scale)) {
      # lm replaces with total variance, we will "scale" by a certain constant
      out$se.fit <- out$se.fit * scale
      df <- df
    } else {
      df <- Inf
    }
  }
  if (interval == "confidence") {
    tstar <- qt(1 - (1 - level) / 2, df = df)
    out$lwr <- out$fit - tstar * se
    out$upr <- out$fit + tstar * se
  }
  # when interval = "prediction", predict.lm instead stores
  # se <- sqrt(se^2 + residual.scale^2) to each fit
  # Not sure if this makes sense for spatial data, so not including it at this point
  # and hence, the prior error for interval = "prediction"
  # as an aside, does it even make sense for lm data? Adding observation variance
  # to the variance of the fixed effects and then creating an interval for
  # the fixed effect itself makes the interval is wider than it should be

  # subset out to include terms objects
  # this is inefficient, and subsetting should be done prior
  # however, given this is only the size of the predictor space, it should
  # not be too burdensome
  if (!is.null(terms)) {
    if (is.numeric(terms) || is.character(terms)) {
      out <- lapply(out, function(x) x[, terms, drop = FALSE])
    } else {
      stop("terms must be character or numeric.", call. = FALSE)
    }
  }
  # return terms
  # remove list structure if only asked for fit
  if (length(out) == 1 && names(out) == "fit") {
    out <- out$fit
  }
  out
}
