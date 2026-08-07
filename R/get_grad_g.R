# the gradient (wrt the free covariance parameters theta) of
# g(theta) = Li' Cov(betahat) Li -- the quantity get_satterthwaite_df()
# needs to delta-method Var(g_hat) from Cov(theta_hat). Cov(betahat) =
# (X'Sigma^-1X)^-1, so this is really differentiating a matrix
# inverse; "numeric" finite-differences get_grad_gi() (which recomputes that
# inverse from scratch at a perturbed theta) directly, while "closed" uses
# the closed-form matrix-derivative identity for (X'Sigma^-1X)^-1.
get_grad_g <- function(Li, method, context, object) {
  UseMethod("get_grad_g", object)
}

#' @noRd
#' @exportS3Method
get_grad_g.splm <- function(Li, method, context, object) {

  if (method == "numeric") {
    # holds every covariance quantity fixed except the free parameters
    # themselves, then lets numDeriv perturb those and re-derive g(theta)
    # from scratch each time (get_grad_gi()) -- works for any covariance
    # type, at the cost of finite-difference approximation error 
    # (especially useful for models with anisotropy)
    obj_grad <- function(cov_val_free) {
      get_grad_gi(cov_val_free, Li, context, object)
    }
    grad_g <- numDeriv::grad(obj_grad, context$cov_val_free)
  } else if (method == "closed") {
    # d/dtheta_k[(X'Sigma^-1X)^-1] = -(X'Sigma^-1X)^-1 X'Sigma^-1 (dSigma/dtheta_k)
    # Sigma^-1 X (X'Sigma^-1X)^-1 (the standard derivative of a matrix
    # inverse, applied to Var(betahat)); contracting both sides on the left
    # and right by Li collapses the two copies of (X'Sigma^-1X)^-1 X'Sigma^-1
    # into a single vector a = Sigma^-1 X Var(betahat) Li, so each entry of
    # the gradient is just the quadratic form a' (dSigma/dtheta_k) a -- exact
    # and far cheaper than differentiating the full p x p matrix inverse
    dSig_list <- get_dSig_dtheta_cov(context, object)
    X <- context$data_object$X_list[[1]]
    a <- chol2inv(chol(covmatrix(object))) %*% X %*% vcov(object) %*% Li
    grad_g <- vapply(dSig_list, function(d_Sigi) {
      as.numeric(crossprod(a, d_Sigi %*% a))
    }, numeric(1))
  }
  names(grad_g) <- context$cov_names_free
  grad_g
}

#' @noRd
#' @exportS3Method
# same intuition as get_grad_g.splm() above -- only the covariance/precision
# machinery underneath (get_grad_gi.spautor()/get_dSig_dtheta_cov.spautor())
# differs to reflect spautor's areal (CAR/SAR) structure
get_grad_g.spautor <- function(Li, method, context, object) {

  if (method == "numeric") {
    obj_grad <- function(cov_val_free) {
      get_grad_gi(cov_val_free, Li, context, object)
    }
    grad_g <- numDeriv::grad(obj_grad, context$cov_val_free)
  } else if (method == "closed") {
    # placeholder: see get_vcov_theta.spautor()'s "closed" branch -- not
    # currently used for car/sar because closed-form derivitaves are not yet
    # derived, keep structure for future updates
    # dSig_list <- get_dSig_dtheta_cov(context, object)
    # X <- context$data_object$X
    # a <- chol2inv(chol(covmatrix(object))) %*% X %*% vcov(object) %*% Li
    # grad_g <- vapply(dSig_list, function(d_Sigi) {
    #   as.numeric(crossprod(a, d_Sigi %*% a))
    # }, numeric(1))
  }
  names(grad_g) <- context$cov_names_free
  grad_g
}

# g(theta) = Li' Cov(betahat) Li itself evaluated at an arbitrary
# (possibly finite-difference-perturbed) free-parameter vector cov_val_free
# -- this is exactly what numDeriv::grad() calls repeatedly at nearby points
# to approximate get_grad_g()'s "numeric" gradient; every fixed/known
# parameter is held at context's original fitted value throughout
# (fill_cov_params() below), so only the free parameters actually move
get_grad_gi <- function(cov_val_free, Li, context, object) {
  UseMethod("get_grad_gi", object)
}

#' @noRd
#' @exportS3Method
get_grad_gi.splm <- function(cov_val_free, Li, context, object) {
  names(cov_val_free) <- context$cov_names_free
  filled <- fill_cov_params(cov_val_free, context)

  if (context$anisotropy) {
    dist_matrix <- as.matrix(build_anis_dist_matrix_list(filled$spcov_params, context$data_object)[[1]])
  } else {
    dist_matrix <- context$dist_matrix
  }
  randcov_Zs <- if (is.null(context$data_object$randcov_list)) NULL else context$data_object$randcov_list[[1]]
  partition_matrix_val <- if (is.null(context$data_object$partition_list)) NULL else context$data_object$partition_list[[1]]

  Sig <- cov_matrix(
    filled$spcov_params, dist_matrix, filled$randcov_params, randcov_Zs,
    partition_matrix_val,
    diagtol = context$data_object$diagtol
  )
  Sig_lowchol <- t(chol(Sig))
  SqrtSigInv_X <- forwardsolve(Sig_lowchol, context$data_object$X_list[[1]])
  vcov_betahat <- chol2inv(chol(forceSymmetric(crossprod(SqrtSigInv_X, SqrtSigInv_X))))
  as.numeric(crossprod(Li, vcov_betahat) %*% Li)
}

#' @noRd
#' @exportS3Method
get_grad_gi.spautor <- function(cov_val_free, Li, context, object) {
  names(cov_val_free) <- context$cov_names_free
  filled <- fill_cov_params(cov_val_free, context)

  # car/sar parameterize the precision matrix directly and sparsely;
  # spautor_cov_matrixInv() (the same helper spautor()'s own likelihood
  # uses, see gloglik_products.car()) already handles M, random effects,
  # partitioning, and reducing the full-W-sized result down to the
  # observed rows via Sherman-Morrison-Woodbury -- reused as-is rather than
  # re-deriving that reduction here
  X <- context$data_object$X
  SigInv <- spautor_cov_matrixInv(
    filled$spcov_params, context$data_object, context$data_object$W,
    filled$randcov_params,
    ldet = FALSE
  )$SigInv
  vcov_betahat <- chol2inv(chol(forceSymmetric(crossprod(X, SigInv %*% X))))
  as.numeric(crossprod(Li, vcov_betahat) %*% Li)
}

# reassembles a full spcov/randcov parameter set by overwriting just the
# free entries of the originally-fitted values with cov_val_free -- known/
# fixed parameters (e.g. a user-supplied spcov_initial(..., known = "de"))
# are left untouched, since only the free ones are ever perturbed by
# numDeriv or differentiated against
fill_cov_params <- function(cov_val_free, context) {
  spcov_params_val <- context$spcov_params
  for (nm in context$spcov_names_free) spcov_params_val[[nm]] <- cov_val_free[[nm]]

  randcov_params_val <- context$randcov_params
  for (nm in context$randcov_names_free) randcov_params_val[[nm]] <- cov_val_free[[nm]]

  list(spcov_params = spcov_params_val, randcov_params = randcov_params_val)
}
