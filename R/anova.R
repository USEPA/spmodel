#' Compute analysis of variance and likelihood ratio tests of fitted model objects
#'
#' @description Compute analysis of variance tables for a fitted model object or
#'   a likelihood ratio test for two fitted model objects.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... An additional fitted model object.
#' @param test A logical value indicating whether p-values from asymptotic Chi-squared
#'   hypothesis tests should be returned. Defaults to \code{TRUE}.
#' @param Terms An optional character or integer vector that specifies terms in the model
#'   used to jointly compute test statistics and p-values (if \code{test = TRUE})
#'   against a null hypothesis of zero. \code{Terms} is only used when a single fitted model
#'   object is passed to the function. If \code{Terms} is a character vector, it
#'   should contain the names of the fixed effect terms. If \code{Terms} is an integer
#'   vector, it should correspond to the order (starting at one) of the names
#'   of the fixed effect terms. The easiest way to obtain the names of
#'   all possible terms is to run \code{tidy(anova(object))$effects} (the
#'   integer representation matches the positions of this vector).
#' @param L An optional numeric matrix or list specifying linear combinations
#'   of the coefficients in the model used to compute test statistics
#'   and p-values (if \code{test = TRUE}) for coefficient constraints corresponding to a null
#'   hypothesis of zero. \code{L} is only used when a single fitted model
#'   object is passed to the function. If \code{L} is a numeric matrix, its rows
#'   indicate coefficient constraints and its columns
#'   represent coefficients. Then a single hypothesis test is conducted
#'   against a null hypothesis of zero.
#'   If \code{L} is a list, each list element is a numeric matrix specified as above.
#'   Then separate hypothesis tests are conducted. The easiest
#'   way to obtain all possible coefficients is to run \code{tidy(object)$term}.
#'
#' @details When one fitted model object is present, \code{anova()}
#'   performs a general linear hypothesis test corresponding to some hypothesis
#'   specified by a matrix of constraints. If \code{Terms} and \code{L} are not specified,
#'   each model term is tested against zero (which correspond to type III or marginal
#'   hypothesis tests from classical ANOVA). If \code{Terms} is specified and \code{L}
#'   is not specified, all terms are tested jointly against zero. When \code{L} is
#'   specified, the linear combinations of terms specified by \code{L} are jointly
#'   tested against zero.
#'
#'   When two fitted model objects are present, one must be a "reduced"
#'   model nested in a "full" model. Then \code{anova()} performs a likelihood ratio test.
#'
#' @return When one fitted model object is present, \code{anova()}
#'   returns a data frame with degrees of
#'   freedom (\code{Df}), test statistics (\code{Chi2}), and p-values
#'   (\code{Pr(>Chi2)} if \code{test = TRUE}) corresponding
#'   to asymptotic Chi-squared hypothesis tests for each model term.
#'
#'   When two fitted model objects are present, \code{anova()} returns a data frame
#'   with the difference in degrees of freedom between the full and reduced model (\code{Df}), a test
#'   statistic (\code{Chi2}), and a p-value corresponding to the likelihood ratio test
#'   (\code{Pr(>Chi2)} if \code{test = TRUE}).
#'
#'   Whether one or two fitted model objects are provided,
#'   \code{tidy()} can be used
#'   to obtain tidy tibbles of the \code{anova(object)} output.
#'
#' @name anova.spmodel
#' @method anova splm
#' @order 1
#' @export
#'
#' @examples
#' # one-model anova
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' anova(spmod)
#' tidy(anova(spmod))
#' # see terms
#' tidy(anova(spmod))$effects
#' tidy(anova(spmod, Terms = c("water", "tarp")))
#' # same as
#' tidy(anova(spmod, Terms = c(2, 3)))
#' # likelihood ratio test
#' lmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "none"
#' )
#' tidy(anova(spmod, lmod))
anova.splm <- function(object, ..., test = TRUE, Terms, L) {
  # see if one or two models
  object2_list <- list(...)

  if (missing(Terms)) Terms <- NULL
  if (missing(L)) L <- NULL

  # one model stuff
  if (length(object2_list) == 0) {
    L <- get_L(L, Terms, object)
    # run the Wald chi-squared test for each hypothesis matrix and stack results
    anova_val <- do.call(rbind, lapply(L, get_marginal_Chi2, object))

    if (!test) {
      anova_val <- anova_val[-which(colnames(anova_val) == "Pr(>Chi2)")]
    }
    anova_val <- structure(anova_val, heading = c("Analysis of Variance Table\n", paste("Response:", deparse(object$formula[[2L]]))))
  }

  # two model stuff
  # likelihood ratio test between a nested pair of models
  else {
    object2 <- object2_list[[1]]
    if (!object$estmethod %in% c("ml", "reml") || !object2$estmethod %in% c("ml", "reml")) {
      stop("LRT only defined for ml or reml", call. = FALSE)
    }

    # reml log-likelihoods are only comparable when the fixed effects are
    # identical (reml profiles out fixed effects, so mixing estmethods or
    # varying fixed effects under reml makes the likelihoods non-comparable)
    if (all(c("ml", "reml") %in% c(object$estmethod, object2$estmethod))) {
      stop("Both fitted model objects must have the same estimation method", call. = FALSE)
    }

    if (
      (object$estmethod %in% c("reml") && object2$estmethod %in% c("reml")) &&
        any(sort(colnames(model.matrix(object))) != sort(colnames(model.matrix(object2))))
    ) {
      stop("The fixed effect coefficients must be the same when performing a likeihood ratio test using the reml estimation method. To perform the likelihood ratio tests for different fixed effect and covariance coefficients simultaneously, refit the models using the ml estimation method.", call. = FALSE)
    }
    # LRT statistic: -2 * (loglik of reduced model - loglik of full model),
    # asymptotically chi-squared under the null that the reduced model holds
    Chi2_stat <- abs(-2 * (logLik(object2) - logLik(object)))

    # df for ml vs reml
    # ml estimates fixed effects + covariance params, reml estimates only
    # covariance params (see AICc.R for the same distinction)
    df1 <- object$npar
    df2 <- object2$npar
    if (object$estmethod == "ml") df1 <- df1 + object$p
    if (object2$estmethod == "ml") df2 <- df2 + object2$p
    df_diff <- abs(df1 - df2)
    p_value <- pchisq(Chi2_stat, df_diff, lower.tail = FALSE)
    # the model with more estimated parameters (npar) is the "full" model;
    # the other is "reduced" -- used only for labeling the output
    if (object2$npar < object$npar) {
      full_name <- deparse(substitute(object)) # replace as.character with deparse
      reduced_name <- as.character(as.list(substitute(list(...)))[-1])
    } else {
      reduced_name <- deparse(substitute(object)) # replace as.character with deparse
      full_name <- as.character(as.list(substitute(list(...)))[-1])
    }
    if (test) {
      anova_val <- data.frame(Df = df_diff, Chi2 = Chi2_stat, p.value = p_value)
      colnames(anova_val) <- c("Df", "Chi2", "Pr(>Chi2)")
    } else {
      anova_val <- data.frame(Df = df_diff, Chi2 = Chi2_stat)
      colnames(anova_val) <- c("Df", "Chi2")
    }
    rownames(anova_val) <- paste(full_name, "vs", reduced_name)
    attr(anova_val, "full") <- full_name
    attr(anova_val, "reduced") <- reduced_name
    anova_val <- structure(anova_val, heading = c("Likelihood Ratio Test\n", paste("Response:", deparse(object$formula[[2L]]))))
  }
  structure(anova_val, class = c(paste("anova", class(object), sep = "."), "data.frame"))
}

#' @rdname anova.spmodel
#' @method anova spautor
#' @order 2
#' @export
anova.spautor <- anova.splm

#' Compute a marginal Wald chi-squared test from a general linear hypothesis matrix
#'
#' @param L A hypothesis matrix (or vector, coerced to a single-row matrix)
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()]
#'
#' @return A data frame with columns \code{Df}, \code{Chi2}, and \code{Pr(>Chi2)}
#'   for the general linear hypothesis test \eqn{L\beta = 0}
#'
#' @noRd
get_marginal_Chi2 <- function(L, object) {
  # make matrix if a numeric vector
  if (!is.matrix(L)) {
    L <- matrix(L, nrow = 1)
  }
  # find the number of rows
  Df <- NROW(L)
  # find product2 of the GLHT
  # inverse of the sampling covariance of L %*% beta_hat, via Cholesky for
  # numerical stability/efficiency rather than a direct matrix inverse
  part2 <- chol2inv(chol(forceSymmetric(L %*% vcov(object) %*% t(L))))
  # find product3 of the GLHT
  part3 <- L %*% coefficients(object)
  # compute the chi-squared statistic
  # Wald statistic (L*beta_hat)' [L*Var(beta_hat)*L']^-1 (L*beta_hat), which
  # is asymptotically chi-squared with Df degrees of freedom under H0: L*beta = 0
  # Chi2/rank(L) is an F(rank(L), Inf) distribution, which equals a scaled chi-squared
  # multiply the F value Chi2/rank(L) by rank(L) yields the original chi-squared
  # with rank(L) df
  Chi2 <- as.numeric(crossprod(part3, part2) %*% part3)
  # find the p-value
  p.value <- pchisq(Chi2, Df, lower.tail = FALSE)
  # put it all in a data frame
  Chi2_df <- data.frame(Df, Chi2, p.value)
  # assign column and row names
  colnames(Chi2_df) <- c("Df", "Chi2", "Pr(>Chi2)")
  rownames(Chi2_df) <- names(L)
  # return the data frame
  Chi2_df
}

#' @rdname anova.spmodel
#' @param x An object from \code{anova(object)}.
#' @method tidy anova.splm
#' @order 5
#' @export
tidy.anova.splm <- function(x, ...) {
  if (!is.null(attr(x, "full")) && !is.null(attr(x, "reduced"))) {
    result <- tibble::tibble(full = attr(x, "full"), reduced = attr(x, "reduced"), df = x$Df, statistic = x$Chi2)
  } else {
    result <- tibble::tibble(effects = rownames(x), df = x$Df, statistic = x$Chi2)
  }
  if ("Pr(>Chi2)" %in% colnames(x)) {
    result$p.value <- x[["Pr(>Chi2)"]]
  }
  result
}

#' @rdname anova.spmodel
#' @method tidy anova.spautor
#' @order 6
#' @export
tidy.anova.spautor <- tidy.anova.splm
