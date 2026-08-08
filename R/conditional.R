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
#'   \code{c("newdata", "beta", "object")}. The default is \code{"newdata"}.
#'   See Details for more.
#' @param type For \code{spglm()} model objects, the scale of the conditional
#'   simulations for \code{newdata}.
#'   When \code{type = "link"}, the predicted means on the
#'   link scale are returned. When \code{type = "response"}, the predicted means
#'   on the response scale are returned. When \code{type = "new"}, a new observation
#'   is simulated from the appropriate response distribution with mean equal to
#'   the mean on the response scale and dispersion equal to the dispersion parameter
#'   from \code{object}. The default is \code{"link"}.
#' @param samples  The number of conditional simulations. The default is
#'   \code{10,000}.
#' @param local An optional logical or list controlling the big data approximation.
#'   If omitted, \code{local} is set
#'   to \code{TRUE} or \code{FALSE} based on the whether the observed
#'   or prediction sample size (the number of
#'   non-missing observations in \code{data} or \code{newdata}) -- if either exceeds 5,000,
#'   \code{local} is set to \code{TRUE}. Otherwise it is set to \code{FALSE}.
#'   If \code{FALSE}, no big data approximation is implemented.
#'   If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{reorder_base}: The data reordering approached to reorder
#'       the observed data prior to subsetting to obtain a base sample. If \code{reorder = "none"}, no reordering
#'       is applied to the observed data. If \code{reorder = "random"}, the observed data order is
#'       randomly reshuffled. If \code{reorder = "grts"}, the observed data order is
#'       randomly generated using the GRTS algorithm for spatially balanced
#'       sampling via \code{spsurvey::grts()}. The default is \code{"grts"}.
#'     \item \code{size_base}: The number of observed data observations used for the base sample.
#'       The default is 3,000. See Details for more.
#'     \item \code{reorder_new}: The data reordering approached to reorder
#'       \code{newdata}. If \code{reorder = "none"}, no reordering
#'       is applied to \code{newdata}. If \code{reorder = "random"}, \code{newdata} is
#'       randomly reshuffled. If \code{reorder = "grts"}, \code{newdata} is
#'       randomly generated using the GRTS algorithm for spatially balanced
#'       sampling via \code{spsurvey::grts()}. The default is \code{"grts"}.
#'     \item \code{kmeans_new}: Whether \code{newdata} observations
#'       they should be assigned to blocks based on k-means clustering
#'       on the coordinates, with clusters of size approximately equal to
#'       \code{size_new}. The default is \code{FALSE} when \code{reorder_new = "none"}
#'       and \code{TRUE} otherwise.
#'     \item \code{size_new}: The (approximate) number of observations used
#'       for each block. The default is 500. See Details for more.
#'       The default is 500.
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(reoder_base = "grts", size_base = 3000, reorder_new = "grts",
#'   size_new = 500, kmeans = TRUE, parallel = FALSE)}.
#
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
#'   fixed effects.
#'
#'   \code{local} Details: The big data approximation works by assigning \code{size_base}
#'   observations to a base sample and then simulating data for the base sample.
#'   The remaining observations are assigned to blocks. For each block, data
#'   are simulated from the conditional distribution given the base sample.
#'   Observations from the same block share conditional covariance while
#'   observations from distinct blocks are assumed conditionally independent
#'   (given the base sample). Parallelization generally further speeds up
#'   computations.
#'
#' @return If \code{output = "newdata"}, an a x b matrix of conditional simulations
#'   for each row in \code{newdata}, where a is the
#'   number of rows in \code{newdata} and b is the number of samples.
#'   If \code{output = "beta"}, an p x b matrix of conditional simulations for each
#'   element in \code{coef(object)}, where p is the
#'   number of fixed effects and b is the number of samples.
#'   If \code{output = "object"}, an n x b matrix of observed data values, where n is the
#'   number of rows in \code{data} and b is the number of samples. This is useful
#'   when the goal is bind together observed data and conditional simulations for
#'   \code{newdata}, as they have the same column dimension.
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
conditional.splm <- function(object, newdata, output = "newdata", samples = 10000, local, ...) {

  if (missing(local)) {
    local <- NULL
  }

  if ("all" %in% output) {
    output <- c("newdata", "beta", "object")
  }
  if (any(!output %in% c("newdata", "beta", "object"))) {
    stop("output must be \"newdata\", \"beta\", \"object\", or \"all\".", call. = FALSE)
  }

  # error if newdata missing from arguments and object
  if (missing(newdata)) {
    if (is.null(object$newdata)) {
      stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
    } else {
      newdata <- object$newdata
    }
  }

  y <- model.response(model.frame(object))
  base_val_y <- matrix(rep(y, times = samples), ncol = samples)
  if (length(output) == 1 && output == "object") {
    return(base_val_y)
  }
  # handle offset
  offset_obdata <- model.offset(model.frame(object))
  if (!is.null(offset_obdata)) {
    y <- y - offset_obdata
  }

  local_list <- get_local_list_conditional(local, object, newdata)

  betahat <- coef(object)
  X <- model.matrix(object)


  cov_betahat_lowchol <- t(chol(vcov(object)))
  new_betahat <- vapply(seq_len(samples), function(x) as.numeric(cov_betahat_lowchol %*% rnorm(length(betahat))), numeric(length(betahat)))
  # beta0 force to matrix
  if (!is.matrix(new_betahat)) {
    new_betahat <- matrix(new_betahat, nrow = 1)
  }
  new_betahat <- sweep(new_betahat, 1, betahat, "+")
  new_fitted <- X %*% new_betahat
  new_resid <- sweep(-1 * new_fitted, 1, y, "+")

  # now simulate beta and add
  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
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
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  if (local_list$method_base != "all") {
    object$obdata <- object$obdata[local_list$index$base, , drop = FALSE]
    base_val <- new_resid[local_list$index$base, , drop = FALSE]
  } else {
    base_val <- new_resid
  }
  if (local_list$method_new != "all") {
    x0 <- lapply(local_list$index$new, function(x) newdata_model[x, , drop = FALSE])
    newdata <- lapply(local_list$index$new, function(x) newdata[x, , drop = FALSE])
  } else {
    x0 <- list(newdata_model)
    newdata <- list(newdata)
  }
  newdata_list <- mapply(x = x0, y = newdata, FUN = function(x, y) list(x0 = x, newdata = y), SIMPLIFY = FALSE)
  spcov_val <- coef(object, type = "spcov")
  if (spcov_val[["de"]] == 0 && is.null(coef(object, type = "randcov"))) {
    cov_lowchol_base <- Matrix::Diagonal(n = object$n, x = sqrt(spcov_val[["ie"]]))
  } else {
    cov_lowchol_base <- t(chol(covmatrix(object)))
  }
  SqrtSigInv_X <- forwardsolve(cov_lowchol_base, X)
  cov_betahat <- vcov(object)
  if (local_list$parallel) {
    cl <- parallel::makeCluster(local_list$ncores)
    new_val <- parLapply(cl, newdata_list, get_conditional_new_from_base_adjust, object, base_val, cov_lowchol_base, samples, SqrtSigInv_X, cov_betahat)
    cl <- parallel::stopCluster(cl)
  } else {
    new_val <- lapply(newdata_list, get_conditional_new_from_base_adjust, object, base_val, cov_lowchol_base, samples, SqrtSigInv_X, cov_betahat)
  }


  new_val <- do.call("rbind", new_val)
  if (local_list$method_new != "all") {
    index_new <- do.call("c", local_list$index$new)
    new_val <- new_val[order(index_new), , drop = FALSE]
  }

  new_val <- newdata_model %*% new_betahat + new_val

  if (!is.null(offset)) {
    new_val <- sweep(new_val, 1, offset, "+")
  }

  val <- list(newdata = new_val, beta = new_betahat, object = base_val_y)
  if (length(output) == 1) {
    return(val[[output]])
  } else {
    return(val[output])
  }
}

#' @rdname conditional
#' @method conditional spglm
#' @export
conditional.spglm <- function(object, newdata, output = "newdata", type = c("link", "response", "new"), samples = 10000, local, newdata_size, ...) {

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

  w <- fitted(object, type = "link")
  y <- object$y
  size <- object$size
  base_val_w <- matrix(rep(w, times = samples), ncol = samples)
  if (length(output) == 1 && output == "object") {
    return(base_val_w)
  }
  # handle offset
  offset_obdata <- model.offset(model.frame(object))
  if (!is.null(offset_obdata)) {
    w <- w - offset_obdata
  }


  local_list <- get_local_list_conditional(local, object, newdata)

  betahat <- coef(object)
  X <- model.matrix(object)

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
  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
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
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  if (local_list$method_base != "all") {
    object$obdata <- object$obdata[local_list$index$base, , drop = FALSE]
    X <- X[local_list$index$base, , drop = FALSE]
    w <- w[local_list$index$base]
    y <- y[local_list$index$base]
    if (!is.null(size)) {
      size <- size[local_list$index$base]
    }
  }

  cov_lowchol_base <- t(chol(covmatrix(object)))
  SigInv <- chol2inv(t(cov_lowchol_base))
  SqrtSigInv_X <- forwardsolve(cov_lowchol_base, X)
  SigInv_X <- backsolve(t(cov_lowchol_base), SqrtSigInv_X)
  cov_betahat <- vcov(object, var_correct = FALSE)
  Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
  D <- get_D(object$family, w, y, size, as.vector(object$coefficients$dispersion))
  cov_lowchol_mH <- t(chol(Matrix::forceSymmetric(-1 * (D - Ptheta)))) # this is actually the inverse of covariance matrix of w
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)

  # cov_lowchol_w <- t(chol(chol2inv(t(cov_lowchol_mH)))) # made this more efficient below but not as correct
  # new_w <- vapply(seq_len(samples), function(x) as.numeric(cov_lowchol_w %*% rnorm(length(w))), numeric(length(w)))

  # increase efficiency by exploiting triangular cholesky relationships to reach
  # above solution but with different rnorm() entries (which are then rearranged)
  cov_lowchol_w <- t(forwardsolve(cov_lowchol_mH, Matrix::Diagonal(length(w))))
  reshuffle <- seq(length(w), 1)
  cov_lowchol_w <- cov_lowchol_w[reshuffle, reshuffle, drop = FALSE]
  new_w <- vapply(seq_len(samples), function(x) as.numeric(cov_lowchol_w %*% rnorm(length(w)))[reshuffle], numeric(length(w)))
  w <- sweep(new_w, 1, w, "+")

  base_val <- w - X %*% new_betahat

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
    new_val <- parLapply(cl, newdata_list, get_conditional_new_from_base_adjust_glm, object, base_val, cov_lowchol_base, samples, SqrtSigInv_X, cov_betahat, SigInv, SigInv_X, wts_beta, cov_lowchol_mH)
    cl <- parallel::stopCluster(cl)
  } else {
    new_val <- lapply(newdata_list, get_conditional_new_from_base_adjust_glm, object, base_val, cov_lowchol_base, samples, SqrtSigInv_X, cov_betahat, SigInv, SigInv_X, wts_beta, cov_lowchol_mH)
  }


  new_val <- do.call("rbind", new_val)
  if (local_list$method_new != "all") {
    index_new <- do.call("c", local_list$index$new)
    new_val <- new_val[order(index_new), , drop = FALSE]
  }

  new_val <- newdata_model %*% new_betahat + new_val

  if (!is.null(offset)) {
    new_val <- sweep(new_val, 1, offset, "+")
  }

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

invlink_conditional <- function(mu_link, type, dispersion, family, newdata_size) {

  # compare to delta method?
  n_newdata <- NROW(mu_link)
  n_sim <- NCOL(mu_link)
  len_sim <- seq(1, n_sim)
  mu <- invlink(mu_link, family, size = 1)
  rm("mu_link")

  if (type == "response") {
    val <- mu
  } else if (type == "new") {

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

  if (type == "response" && family == "binomial") {
    val <- sweep(val, 1, newdata_size, "*")
  }

  val
}
