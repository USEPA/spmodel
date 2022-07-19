#' Compute analysis of variance and likelihood ratio tests of fitted model objects
#'
#' @description Compute analysis of variance tables for a fitted model object or
#'   a likelihood ratio test for two fitted model objects.
#'
#' @param object A fitted model object from [splm()] or [spautor()]
#' @param ... An additional fitted model object from [splm()] or [spautor()]
#'   (for \code{anova()}).
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
#'
#' @method anova spmod
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
anova.spmod <- function(object, ..., test = TRUE, Terms, L) {


  # see if one or two models
  object2_list <- list(...)

  # one model stuff
  if (length(object2_list) == 0) {
    if (missing(L)) {
      assign_indices <- attr(model.matrix(object), "assign") + 1
      # attr(model.matrix(object), "assign") if centering at zero
      if (missing(Terms)) {
        assign_index <- unique(assign_indices)
        L <- lapply(assign_index, get_L_list, assign_indices)
        label <- labels(object)
        if (attr(terms(object), "intercept") == 1) {
          label <- c("(Intercept)", label)
        }
        names(L) <- label
      } else {
        if (is.character(Terms)) {
          Terms <- which(c("(Intercept)", labels(object)) %in% Terms) # - 1 if centering at zero
        }
        L <- list(do.call(rbind, lapply(Terms, get_L_list, assign_indices)))
        label <- c("(Intercept)", labels(object))
        label <- label[Terms] # label[Terms + 1] if centering at zero
        names(L) <- paste(label, collapse = ", ")
      }
    } else {
      if (!is.list(L)) {
        L <- list(L)
      }
      names(L) <- paste("contrast", seq_along(L), sep = "")
    }
    anova_val <- do.call(rbind, lapply(L, get_marginal_Chi2, object))

    if (!test) {
      anova_val <- anova_val[-which(colnames(anova_val) == "Pr(>Chi2)")]
    }
    anova_val <- structure(anova_val, heading = c("Analysis of Variance Table\n", paste("Response:", deparse(object$formula[[2L]]))))
  }

  # two model stuff
  else {
    object2 <- object2_list[[1]]
    if (!object$estmethod %in% c("ml", "reml") || !object2$estmethod %in% c("ml", "reml")) {
      stop("LRT only defined for ml or reml", call. = FALSE)
    }

    if (all(c("ml", "reml") %in% c(object$estmethod, object2$estmethod))) {
      stop("Both fitted model objects must have the same estimation method", call. = FALSE)
    }

    if (
      (object$estmethod %in% c("reml") && object2$estmethod %in% c("reml")) &&
        any(sort(colnames(model.matrix(object))) != sort(colnames(model.matrix(object2))))
    ) {
      stop("The fixed effect coefficients must be the same when performing a likeihood ratio test using the reml estimation method. To perform the likelihood ratio tests for different fixed effect and covariance coefficients simultaneously, refit the models using the ml estimation method.", call. = FALSE)
    }
    Chi2_stat <- abs(-2 * (logLik(object2) - logLik(object)))
    df_diff <- abs(object2$npar - object$npar)
    p_value <- pchisq(Chi2_stat, df_diff, lower.tail = FALSE)
    if (object2$npar < object$npar) {
      full_name <- as.character(substitute(object))
      reduced_name <- as.character(as.list(substitute(list(...)))[-1])
    } else {
      reduced_name <- as.character(substitute(object))
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

    # if (object2$npar < object$npar) {
    #   red <- "1"
    #   full <- "2"
    # } else {
    #   red <- "2"
    #   full <- "1"
    # }
    # rownames(anova_val) <- paste(red, "vs", full)
    anova_val <- structure(anova_val, heading = c("Likelihood Ratio Test\n", paste("Response:", deparse(object$formula[[2L]]))))
  }
  structure(anova_val, class = c(paste("anova", class(object), sep = "."), "data.frame"))
}

get_marginal_Chi2 <- function(L, object) {
  # make matrix if a numeric vector
  if (!is.matrix(L)) {
    L <- matrix(L, nrow = 1)
  }
  # find the number of rows
  Df <- NROW(L)
  # find product2 of the GLHT
  part2 <- chol2inv(chol(forceSymmetric(L %*% vcov(object) %*% t(L))))
  # find product3 of the GLHT
  part3 <- L %*% coefficients(object)
  # compute the chi-squared statistic
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

#' @rdname anova.spmod
#' @param x An object from \code{anova(object)}.
#'
#' @method tidy anova.spmod
#' @export
tidy.anova.spmod <- function(x, ...) {
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

# anova deprecated (old method of smw and sums of squares)

# This was the old method of impelementing ANOVA based on sums of squares that was
# deprecated in favor of the general linear hypothesis test, which more easily
# accommodates joint tests and contrasts.


# get_anova <- function(formula, X, y, Sig_lowchol, ss_f) {
#   assign_index <- unique(attr(X, "assign"))
#   terms_val <- terms(formula)
#   label <- labels(terms_val)
#   if (attr(terms_val, "intercept")) {
#     label <- c("(Intercept)", label)
#   }
#   if (length(label) > 1) {
#     anova_val <- do.call(rbind, lapply(assign_index, get_marginal_X2, X, y, Sig_lowchol, ss_f))
#     colnames(anova_val) <- c("Df", "X2", "Pr(>X2)")
#     rownames(anova_val) <- label
#   } else {
#     anova_val <- NULL
#   }
#   anova_val
# }
#
# # alternative assuming df = Inf
# get_marginal_X2 <- function(assign_index, X, y, Sig_lowchol, ss_f) {
#   newX <- X[, attr(X, "assign") != assign_index, drop = FALSE]
#   sqrtSigInv_newX <- forwardsolve(Sig_lowchol, newX)
#   sqrtSigInv_y <- forwardsolve(Sig_lowchol, y)
#   new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(sqrtSigInv_newX, sqrtSigInv_newX))))
#   new_betahat <- new_covbetahat %*% crossprod(sqrtSigInv_newX, sqrtSigInv_y)
#   sqrtSigInv_r <- sqrtSigInv_y - sqrtSigInv_newX %*% new_betahat
#   ss_r <- crossprod(sqrtSigInv_r, sqrtSigInv_r)
#   ss_diff <- ss_r - ss_f
#   df_diff <- sum(attr(X, "assign") == assign_index)
#   X2 <- as.vector(ss_r - ss_f)
#   data.frame(Df = df_diff, X2 = X2, p.value = pchisq(X2, df_diff, lower.tail = FALSE))
# }

# sum of squares could just be
# betahat where relevant columns of x are held out and SigmaHat is still used
# to compute new betas
# then sums of square error are just (y - x newbeta)t SigmaHat (y - xnewbeta)


# # not assuming df = Inf
# get_marginal_f <- function(assign_index, model1) {
#   newX <- model1$model$X[, attr(model1$model$X, "assign") != assign_index, drop = FALSE]
#   sqrtSigInv_newX <- forwardsolve(model1$Sig_lowchol, newX)
#   sqrtSigInv_y <- forwardsolve(model1$Sig_lowchol, model1$model$y)
#   new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(sqrtSigInv_newX, sqrtSigInv_newX))))
#   new_betahat <- new_covbetahat %*% crossprod(sqrtSigInv_newX, sqrtSigInv_y)
#   sqrtSigInv_r <- sqrtSigInv_y - sqrtSigInv_newX %*% new_betahat
#   ss_r <- crossprod(sqrtSigInv_r, sqrtSigInv_r)
#   ss_f <- deviance(model1)
#   ss_diff <- ss_r - ss_f
#   df_diff <- sum(attr(model1$model$X, "assign") == assign_index)
#   df_f <- model1$n - model1$npar
#   ms_diff <- ss_diff / df_diff
#   ms_f <- ss_f / df_f
#   f <- as.vector(ms_diff / ms_f)
#   data.frame(f = f, numDF = df_diff, denDF = df_f, p.value = pf(f, df_diff, df_f, lower.tail = FALSE))
# }

# assuming df = Inf
# get_marginal_f <- function(assign_index, model1) {
#   newX <- model1$model$X[, attr(model1$model$X, "assign") != assign_index, drop = FALSE]
#   sqrtSigInv_newX <- forwardsolve(model1$Sig_lowchol, newX)
#   sqrtSigInv_y <- forwardsolve(model1$Sig_lowchol, model1$model$y)
#   new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(sqrtSigInv_newX, sqrtSigInv_newX))))
#   new_betahat <- new_covbetahat %*% crossprod(sqrtSigInv_newX, sqrtSigInv_y)
#   sqrtSigInv_r <- sqrtSigInv_y - sqrtSigInv_newX %*% new_betahat
#   ss_r <- crossprod(sqrtSigInv_r, sqrtSigInv_r)
#   ss_f <- deviance(model1)
#   ss_diff <- ss_r - ss_f
#   df_diff <- sum(attr(model1$model$X, "assign") == assign_index)
#   df_f <- Inf
#   ms_diff <- ss_diff / df_diff
#   ms_f <- 1
#   f <- as.vector(ms_diff / ms_f)
#   data.frame(Chisq = f, DF = df_diff, p.value = pf(f, df_diff, df_f, lower.tail = FALSE))
# }
