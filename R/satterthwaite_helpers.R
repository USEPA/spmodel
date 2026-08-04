get_vcov_theta <- function(method, context, object) {

  if (method == "numeric") {
    if (context$anisotropy) {
      obj_eta <- function(eta) {
        gloglik_anis(
          eta, context$spcov_orig2optim, context$data_object, context$estmethod,
          spcov_profiled = FALSE, context$randcov_orig2optim,
          randcov_profiled = FALSE
        )
      }
    } else {
      obj_eta <- function(eta) {
        gloglik(
          eta, context$spcov_orig2optim, context$data_object, context$estmethod,
          list(context$dist_matrix),
          spcov_profiled = FALSE, context$randcov_orig2optim,
          randcov_profiled = FALSE
        )
      }
    }


    H_eta <- numDeriv::hessian(obj_eta, context$eta)

    H_eta <- forceSymmetric(H_eta)
    vcov_eta <- tryCatch(chol2inv(chol(H_eta/2)), error = function(e) NULL)
    if (is.null(vcov_eta)) {
      warning("The covariance matrix of the covariance parameters is not numerically positive definite.", call. = FALSE)
      return(NULL)
    }

    # theta(eta): unpack the optim-scale vector back to the original scale and
    # keep only the free entries (the matrix is technically diagonal but generally computed
    # via a full numerical Jacobian, robust to any covariance type's transform)
    delta_transform_eta <- function(eta) {
      unpacked <- unpack_optim2orig(
        context$spcov_orig2optim, context$randcov_orig2optim, eta,
        FALSE, FALSE, context$data_object
      )
      full <- c(as.numeric(unpacked$spcov_params_val), as.numeric(unpacked$randcov_params_val))
      names(full) <- c(names(unpacked$spcov_params_val), names(unpacked$randcov_params_val))
      full[context$cov_names_free]
    }

    J <- numDeriv::jacobian(delta_transform_eta, context$eta)
    vcov_theta <- J %*% base::tcrossprod(vcov_eta, J)
    dimnames(vcov_theta) <- list(context$cov_names_free, context$cov_names_free)



  } else if (method == "closed") {
    dSig_list <- get_dSig_dtheta_cov(context)
    Sig <- covmatrix(object)
    SigInv <- chol2inv(chol(Sig))
    if (object$estmethod == "reml") {
      SigInv_X <- SigInv %*% context$X
      Ptheta <- SigInv - SigInv_X %*% tcrossprod(vcov(object), SigInv_X)
      P_dSig_list <- lapply(dSig_list, function(dSig_i) Ptheta %*% dSig_i)
    } else if (object$estmethod == "ml") {
      P_dSig_list <- lapply(dSig_list, function(dSig_i) SigInv %*% dSig_i)
    }

    k <- length(dSig_list)
    Imat_expected <- matrix(0, k, k)
    colnames(Imat_expected) <- names(dSig_list)
    rownames(Imat_expected) <- names(dSig_list)

    for (i in seq_len(k)) {
      for (j in seq(i, k)) {
        # tr(A %*% B) = sum(A * t(B))
        Imat_expected[i, j] <- 0.5 * sum(P_dSig_list[[i]] * t(P_dSig_list[[j]]))
        if (i != j) Imat_expected[j, i] <- Imat_expected[i, j]
      }
    }

    vcov_theta <- tryCatch(chol2inv(chol(Imat_expected)), error = function(e) NULL)
    if (is.null(vcov_theta)) {
      warning("The covariance matrix of the covariance parameters is not numerically positive definite.", call. = FALSE)
      return(NULL)
    }
    dimnames(vcov_theta) <- dimnames(Imat_expected)
  }

  vcov_theta
}

validate_satterthwaite_scope <- function(object) {
  if (!object$estmethod %in% c("reml", "ml")) {
    stop("Satterthwaite df are only defined for estmethod \"reml\" or \"ml\".", call. = FALSE)
  }
  if (!is.null(object$local_index) || object$n >= 500) {
    stop("Satterthwaite df can only be used for models with n <= 500 fit without local.", call. = FALSE)
  }
  invisible(TRUE)
}
