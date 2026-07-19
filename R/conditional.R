conditional <- function(object, ...) {
  UseMethod("conditional", object)
}

conditional.splm <- function(object, newdata, type = "newdata", samples = 1, local) {

  if (missing(local)) {
    local <- NULL
  }

  if ("all" %in% type) {
    type <- c("newdata", "beta", "object")
  }
  if (any(!type %in% c("newdata", "beta", "object"))) {
    stop("type must be \"newdata\", \"beta\", \"object\", or \"all\".", call. = FALSE)
  }

  y <- model.response(model.frame(object))
  base_val_y <- matrix(rep(y, times = samples), ncol = samples)
  if (length(type) == 1 && type == "object") {
    return(base_val_y)
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
  cov_lowchol_base <- t(chol(covmatrix(object)))
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

  val <- list(newdata = new_val, beta = new_betahat, object = base_val_y)
  if (length(type) == 1) {
    return(val[[type]])
  } else {
    return(val[type])
  }
}

conditional.spglm <- function(object, newdata, type = "newdata", samples = 1, local) {

  if (missing(local)) {
    local <- NULL
  }

  if ("all" %in% type) {
    type <- c("newdata", "beta", "object")
  }
  if (any(!type %in% c("newdata", "beta", "object"))) {
    stop("type must be \"newdata\", \"beta\", \"object\", or \"all\".", call. = FALSE)
  }

  w <- fitted(object, type = "link")
  y <- object$y
  size <- object$size
  base_val_w <- matrix(rep(w, times = samples), ncol = samples)
  if (length(type) == 1 && type == "object") {
    return(base_val_w)
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
  new_resid <- sweep(-1 * new_fitted, 1, w, "+")

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
    w <- w[local_list$index$base]
    y <- y[local_list$index$base]
    if (!is.null(size)) {
      size <- size[local_list$index$base]
    }
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
  cov_lowchol_base <- t(chol(covmatrix(object)))
  SigInv <- chol2inv(t(cov_lowchol_base))
  SqrtSigInv_X <- forwardsolve(cov_lowchol_base, X)
  SigInv_X <- backsolve(t(cov_lowchol_base), SqrtSigInv_X)
  cov_betahat <- vcov(object, var_correct = FALSE)
  Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
  D <- get_D(object$family, w, y, size, as.vector(object$coefficients$dispersion))
  cov_lowchol_mH <- t(chol(Matrix::forceSymmetric(-1 * (D - Ptheta))))
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
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

  val <- list(newdata = new_val, beta = new_betahat, object = base_val_w)
  if (length(type) == 1) {
    return(val[[type]])
  } else {
    return(val[type])
  }
}
