#' Recover data for \code{emmeans} support
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()]
#' @param frame The model frame
#' @param ... Additional arguments passed to \code{emmeans::recover_data()}
#'
#' @return The recovered data, for use by the \code{emmeans} package
#'
#' @details Registered dynamically for the \code{emmeans} generic in
#'   \code{.onLoad()} rather than exported directly, since \code{emmeans}
#'   is only in Suggests.
#'
#' @noRd
recover_data.splm <- function(object, frame = model.frame(object), ...) {
  # check to see if emmeans installed
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Install the emmeans package before using", call. = FALSE)
  }
  # recover data (using emmeans code)
  fcall <- object$call
  # recognize that lm objects have a $model element that is model.frame(object)
  # delete.response() strips the response from the terms so emmeans can build
  # a reference grid of predictors without needing the (possibly absent) outcome
  emmeans::recover_data(fcall, delete.response(terms(object)), frame = frame, na.action = NULL, ...)
}

recover_data.spautor <- recover_data.splm

recover_data.spglm <- recover_data.splm

recover_data.spgautor <- recover_data.splm


#' Build the \code{emmeans} basis for a fitted model
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()]
#' @param trms Model terms
#' @param xlev Factor levels
#' @param grid A reference grid
#' @param ... Additional arguments passed to \code{emmeans} helpers
#'
#' @return A list describing the linear model basis (design matrix, coefficients,
#'   covariance matrix, etc.), for use by the \code{emmeans} package
#'
#' @details Registered dynamically for the \code{emmeans} generic in
#'   \code{.onLoad()} rather than exported directly, since \code{emmeans}
#'   is only in Suggests.
#'
#' @noRd
emm_basis.splm <- function(object, trms, xlev, grid, ...) {
  # check to see if emmeans installed
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Install the emmeans package before using", call. = FALSE)
  }
  bhat <- coef(object)
  nm <- if (is.null(names(bhat))) {
    row.names(bhat)
  } else {
    names(bhat)
  }
  m <- suppressWarnings(model.frame(trms, grid,
    na.action = na.pass,
    xlev = xlev
  ))
  X <- model.matrix(trms, m, contrasts.arg = object$contrasts)
  assign <- attr(X, "assign")
  # reorder/select columns of X to match the order of the fitted coefficient
  # names so X %*% bhat lines up correctly
  X <- X[, nm, drop = FALSE]
  bhat <- as.numeric(bhat)
  V <- emmeans::.my.vcov(object, ...)
  nbasis <- estimability::all.estble # returns a 1x1 NA which says all functions estimable
  misc <- list()
  # splm/spautor fits use Satterthwaite denominator df (matching object$ddf,
  # see get_emmeans_dffun()) when available; spglm/spgautor (and splm/spautor
  # fits without ddf) fall back to the original asymptotic (Inf) df
  if (inherits(object, c("splm", "spautor"))) {
    dfspec <- get_emmeans_dffun(object)
  } else {
    dfspec <- list(dffun = function(k, dfargs) Inf, dfargs = list(), mesg = "asymptotic")
  }
  dffun <- dfspec$dffun
  dfargs <- dfspec$dfargs
  attr(dffun, "mesg") <- dfspec$mesg
  mm <- model.matrix(object)
  mm <- emmeans::.cmpMM(mm, assign = attr(mm, "assign"))
  if (inherits(object, c("spglm", "spgautor"))) {
    # spmodel doesn't store the link name directly on glm objects, so look it
    # up from the family (binomial/beta use logit, everything else uses log)
    famdat <- data.frame(
      family = c("poisson", "nbinomial", "binomial", "beta", "Gamma", "inverse.gaussian")
    )
    famdat$link <- ifelse(famdat$family %in% c("binomial", "beta"), "logit", "log")
    fam <- famdat[match(object$family, famdat$family), ]
    fam <- list(family = fam$family, link = fam$link)
    misc <- emmeans::.std.link.labels(fam, misc)
  }
  list(
    X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun,
    dfargs = dfargs, misc = misc, model.matrix = mm
  )
}

emm_basis.spautor <- emm_basis.splm

emm_basis.spglm <- emm_basis.splm

emm_basis.spgautor <- emm_basis.splm
