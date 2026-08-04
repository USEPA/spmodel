get_grad_g <- function(Li, method, context, object) {

  if (method == "numeric") {
    obj_grad <- function(cov_val_free) {
      get_grad_gi(cov_val_free, Li, context)
    }
    grad_g <- numDeriv::grad(obj_grad, context$cov_val_free)
  } else if (method == "closed") {
    dSig_list <- get_dSig_dtheta_cov(context)
    a <- chol2inv(chol(covmatrix(object))) %*% context$X %*% vcov(object) %*% Li
    grad_g <- vapply(dSig_list, function(d_Sigi) {
      as.numeric(crossprod(a, d_Sigi %*% a))
    }, numeric(1))
  }
  names(grad_g) <- context$cov_names_free
  grad_g
}

get_grad_gi <- function(cov_val_free, Li, context) {
  names(cov_val_free) <- context$cov_names_free
  filled <- fill_cov_params(cov_val_free, context)
  if (context$anisotropy) {
    dist_matrix <- as.matrix(build_anis_dist_matrix_list(filled$spcov_params, context$data_object)[[1]])
  } else {
    dist_matrix <- context$dist_matrix
  }
  Sig <- cov_matrix(
    filled$spcov_params, dist_matrix, filled$randcov_params, context$randcov_Zs,
    context$partition_matrix,
    diagtol = context$diagtol
  )
  Sig_lowchol <- t(chol(Sig))
  SqrtSigInv_X <- forwardsolve(Sig_lowchol, context$X)
  vcov_betahat <- chol2inv(chol(forceSymmetric(crossprod(SqrtSigInv_X, SqrtSigInv_X))))
  as.numeric(crossprod(Li, vcov_betahat) %*% Li)
}

fill_cov_params <- function(cov_val_free, context) {
  spcov_params_val <- context$spcov_params
  for (nm in context$spcov_names_free) spcov_params_val[[nm]] <- cov_val_free[[nm]]

  randcov_params_val <- context$randcov_params
  for (nm in context$randcov_names_free) randcov_params_val[[nm]] <- cov_val_free[[nm]]

  list(spcov_params = spcov_params_val, randcov_params = randcov_params_val)
}
