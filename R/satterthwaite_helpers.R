get_vcov_theta <- function(method, context, object) {
  UseMethod("get_vcov_theta", object)
}

#' @noRd
#' @exportS3Method
# Cov(thetahat), the sampling covariance matrix of the estimated free
# covariance parameters used by get_satterthwaite_df()
# There are two ways to compute it, selected by `method`:
#  - "numeric": standard MLE asymptotics, Cov(etahat) ~= [Fisher info]^-1,
#    approximated by numerically differentiating the fitting log-likelihood
#    itself at the fitted optimizer-scale value (eta), then delta-method
#    mapped back to the original (theta) parameterization
#  - "closed": the expected (Fisher) information matrix computed in closed
#    form from the known derivatives of Sigma with respect to theta -- exact
#    (no finite-difference error) but only available for covariance types
#    with a dSig_dtheta_spcov.<type>() implementation we have derived
get_vcov_theta.splm <- function(method, context, object) {

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
      # splm's dist_matrix_list is a list with a single element
      obj_eta <- function(eta) {
        gloglik(
          eta, context$spcov_orig2optim, context$data_object, context$estmethod,
          list(context$dist_matrix),
          spcov_profiled = FALSE, context$randcov_orig2optim,
          randcov_profiled = FALSE
        )
      }
    }

    # obj_eta() returns -2*log-likelihood (see gloglik()/gloglik_anis()), so
    # its Hessian is -2 * Hessian(log-lik); dividing by 2 converts that into
    # the observed Fisher information (the negative log-lik Hessian), whose
    # inverse is the usual MLE asymptotic covariance estimate
    #
    # numDeriv perturbs the covariance parameters at many nearby points to
    # build this Hessian via finite differences, and those perturbations
    # routinely wander into numerically awkward (but perfectly legitimate)
    # regions -- e.g. matern's smoothness parameter pushing base R's
    # besselK() to warn "value out of range in 'bessel_k'", or a near-zero
    # variance component briefly going slightly negative. These warnings are
    # incidental to the finite-differencing process itself, not a sign that
    # the resulting Cov(theta_hat) is wrong (spmodel already checks that
    # separately via the positive-definiteness check just below), and a
    # single splm()/spautor() fit can otherwise emit thousands of them.
    # Genuine failures still surface: this suppression only affects
    # warning(), not error(), and stop()'d likelihood evaluations propagate
    # normally.
    H_eta <- suppressWarnings(numDeriv::hessian(obj_eta, context$eta))

    H_eta <- forceSymmetric(H_eta)
    vcov_eta <- tryCatch(chol2inv(chol(H_eta / 2)), error = function(e) NULL)
    if (is.null(vcov_eta)) {
      warning("The covariance matrix of the covariance parameters is not numerically positive definite.", call. = FALSE)
      return(NULL)
    }

    # theta(eta): unpack the optim-scale vector back to the original scale and
    # keep only the free entries (the matrix is technically diagonal but generally computed
    # via a full numerical Jacobian because it accommodates any transformation (e.g., exp vs expit))
    delta_transform_eta <- function(eta) {
      unpacked <- unpack_optim2orig(
        context$spcov_orig2optim, context$randcov_orig2optim, eta,
        FALSE, FALSE, context$data_object
      )
      full <- c(as.numeric(unpacked$spcov_params_val), as.numeric(unpacked$randcov_params_val))
      names(full) <- c(names(unpacked$spcov_params_val), names(unpacked$randcov_params_val))
      full[context$cov_names_free]
    }

    # delta method: Cov(theta_hat) ~= J Cov(eta_hat) J', where J is the
    # Jacobian of the (nonlinear, e.g. log/logit) optim-scale-to-original-scale
    # transform evaluated at the fitted value -- maps the Hessian-based
    # covariance above off of the unconstrained optimizer scale and onto the
    # original, interpretable covariance parameter scale used everywhere else
    #
    # suppressWarnings(): see the H_eta call above -- same finite-difference
    # perturbation, same incidental warnings
    J <- suppressWarnings(numDeriv::jacobian(delta_transform_eta, context$eta))
    vcov_theta <- J %*% base::tcrossprod(vcov_eta, J)
    dimnames(vcov_theta) <- list(context$cov_names_free, context$cov_names_free)
  } else if (method == "closed") {
    # expected (Fisher) information for a Gaussian likelihood in terms of the
    # score's covariance: I_jk = 0.5 * tr(P dSig/dtheta_j P dSig/dtheta_k),
    # where P = Sigma^-1 for ML, and REML's P further accounts for the 
    # marginalization of beta
    dSig_list <- get_dSig_dtheta_cov(context, object)
    Sig <- covmatrix(object)
    SigInv <- chol2inv(chol(Sig))
    if (object$estmethod == "reml") {
      SigInv_X <- SigInv %*% context$data_object$X_list[[1]]
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

    # Cov(thetahat) ~= [expected information]^-1, the closed-form analog of
    # the numeric branch's inverse-Hessian -- no delta method needed here
    # since Imat_expected is already computed directly on the theta scale
    vcov_theta <- tryCatch(chol2inv(chol(Imat_expected)), error = function(e) NULL)
    if (is.null(vcov_theta)) {
      warning("The covariance matrix of the covariance parameters is not numerically positive definite.", call. = FALSE)
      return(NULL)
    }
    dimnames(vcov_theta) <- dimnames(Imat_expected)
  }

  vcov_theta
}

#' @noRd
#' @exportS3Method
# same intuition as get_vcov_theta.splm() above (Hessian-of-the-likelihood ->
# Fisher information -> delta method onto the theta scale for "numeric";
# closed-form expected information for "closed") -- spautor only differs in
# how the areal (CAR/SAR) likelihood/W matrix get passed through
get_vcov_theta.spautor <- function(method, context, object) {

  if (method == "numeric") {
    # spautor is never anisotropic, and its dist_matrix_list convention is
    # the raw W matrix itself, not wrapped in a list (see
    # get_data_object_spautor()/get_cov_gradients_context.spautor())
    obj_eta <- function(eta) {
      gloglik(
        eta, context$spcov_orig2optim, context$data_object, context$estmethod,
        context$data_object$W,
        spcov_profiled = FALSE, context$randcov_orig2optim,
        randcov_profiled = FALSE
      )
    }

    # suppressWarnings(): numDeriv's finite-difference perturbations of the
    # covariance parameters routinely wander into numerically awkward (but
    # legitimate) regions -- e.g. besselK() range warnings for matern-family
    # covariances -- that are incidental to the differencing process itself,
    # not evidence the resulting Cov(theta_hat) is wrong (checked separately
    # via the positive-definiteness check just below); a single fit can
    # otherwise emit thousands of these. error()s still propagate normally.
    H_eta <- suppressWarnings(numDeriv::hessian(obj_eta, context$eta))

    H_eta <- forceSymmetric(H_eta)
    vcov_eta <- tryCatch(chol2inv(chol(H_eta / 2)), error = function(e) NULL)
    if (is.null(vcov_eta)) {
      warning("The covariance matrix of the covariance parameters is not numerically positive definite.", call. = FALSE)
      return(NULL)
    }

    # theta(eta): unpack the optim-scale vector back to the original scale and
    # keep only the free entries (the matrix is technically diagonal but generally computed
    # via a full numerical Jacobian because it accommodates any transformation (e.g., exp vs expit))
    delta_transform_eta <- function(eta) {
      unpacked <- unpack_optim2orig(
        context$spcov_orig2optim, context$randcov_orig2optim, eta,
        FALSE, FALSE, context$data_object
      )
      full <- c(as.numeric(unpacked$spcov_params_val), as.numeric(unpacked$randcov_params_val))
      names(full) <- c(names(unpacked$spcov_params_val), names(unpacked$randcov_params_val))
      full[context$cov_names_free]
    }

    # suppressWarnings(): see the H_eta call above -- same finite-difference
    # perturbation, same incidental warnings
    J <- suppressWarnings(numDeriv::jacobian(delta_transform_eta, context$eta))
    vcov_theta <- J %*% base::tcrossprod(vcov_eta, J)
    dimnames(vcov_theta) <- list(context$cov_names_free, context$cov_names_free)
  } else if (method == "closed") {
    # placeholder: no dSig_dtheta_spcov.car()/.sar() exist yet, so this
    # errors informatively (via get_dSig_dtheta_cov()'s default dispatch)
    # rather than doing anything -- get_satterthwaite_method() already
    # forces method = "numeric" for car/sar, so this branch is currently used
    # but is kept so the two vcov_theta methods stay
    # structurally similar and closed-form support can be added later without
    # major structural changes
    dSig_list <- get_dSig_dtheta_cov(context, object)
    Sig <- covmatrix(object)
    SigInv <- chol2inv(chol(Sig))
    if (object$estmethod == "reml") {
      SigInv_X <- SigInv %*% context$data_object$X
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

# Validate elements of the object and method to make sure the parent function can continue:
#  - estmethod reml/ml: the expected-information and log-likelihood-Hessian
#    derivations above assume betahat/theta_hat come from (RE)ML, so other
#    estmethods (e.g. sv-wls) have no likelihood for these to differentiate
#  - no local (big-data approximation): local fitting partitions the data and
#    approximates Sigma piecewise, so there is no single well-defined
#    likelihood/covariance matrix left to differentiate globally
#  - n >= 500: not a correctness issue, just a performance one -- the
#    numerical Hessian/Jacobian calls are O(n^3)-per-evaluation and this can
#    get slow, so users are warned
validate_satterthwaite_scope <- function(object, method) {

    if (method == "numeric" && !requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Install the numDeriv package before using satterthwaite(method = \"numeric\").", call. = FALSE)
  }
  if (!method %in% c("closed", "numeric")) stop("method must be \"closed\" or \"numeric\".", call. = FALSE)

  if (!object$estmethod %in% c("reml", "ml")) {
    stop("Satterthwaite df are only defined for estmethod \"reml\" or \"ml\".", call. = FALSE)
  }
  if (!is.null(object$local_index)) {
    stop("Satterthwaite df can only be used for models fit without local.", call. = FALSE)
  }
  if (object$n >= 500) {
    warning("For sample size n >= 500, Satterthwaite ddf may result in exceedingly long computational times. Consider using the asymptotic ddf instead, which should be similar given the sample size.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Determine a \code{ddf} argument's method
#'
#' Determine whether to use satterthwaite or asymptotic degrees of freedom
#'
#' @param ddf The \code{ddf} argument as passed by the caller (already
#'   resolved from \code{missing()} to \code{NULL})
#' @param n The fitted model's sample size
#'
#' @return Either \code{"asymptotic"} or \code{"satterthwaite"}
#'
#' @noRd
determine_ddf <- function(ddf, n) {
  if (is.null(ddf)) {
    ddf <- if (n <= 500) "satterthwaite" else "asymptotic"
  }
  if (! ddf %in% c("asymptotic", "satterthwaite")) {
    stop("ddf must be \"asymptotic\" or \"satterthwaite\".", call. = FALSE)
  }
  ddf
}

#' Compute (or skip) denominator degrees of freedom for a fitted model
#'
#' Implements \code{splm()}/\code{spautor()}'s \code{ddf} argument (see
#' \code{determine_ddf()} for how it resolves the sample size to
#' \code{"asymptotic"} or \code{"satterthwaite"}): \code{"asymptotic"} returns
#' every element \code{NULL} (unchanged from before this argument existed);
#' \code{"satterthwaite"} attempts \code{satterthwaite_core()} on the
#' already-fitted object and returns everything it produces, or every
#' element \code{NULL} if the computation errors for any reason -- a failure
#' here must never prevent \code{splm()}/\code{spautor()} from returning a
#' fitted model. \code{vcov_cov}/\code{vcov_spcov}/\code{vcov_randcov} (the
#' estimated covariance matrix of the free covariance parameters, as a whole
#' and split into its spatial/random-effect blocks) are byproducts of
#' computing \code{ddf} already needed here, so all four are returned
#' together rather than computing \code{vcov_theta} a second time via
#' \code{satterthwaite()} -- see \code{vcov.spmodel()}, which exposes them as
#' \code{vcov(object, type = "cov"/"spcov"/"randcov")}.
#'
#' @param object The fitted model object.
#' @param ddf The \code{ddf} method.
#'
#' @return A list with elements \code{ddf} (a named numeric vector of
#'   denominator df, see \code{satterthwaite()}, or \code{NULL}),
#'   \code{vcov_cov} (the corresponding covariance-parameter covariance
#'   matrix, or \code{NULL}), and \code{vcov_spcov}/\code{vcov_randcov}
#'   (\code{vcov_cov}'s spatial-only/random-effect-only blocks, or
#'   \code{NULL})
#'
#' @noRd
get_fit_ddf <- function(object, ddf) {
  none <- list(ddf = NULL, vcov_cov = NULL, vcov_spcov = NULL, vcov_randcov = NULL)

  ddf <- determine_ddf(ddf, object$n)
  if (ddf == "asymptotic") {
    return(none)
  }

  out <- tryCatch(satterthwaite_core(object, NULL), error = function(e) NULL)
  if (is.null(out)) {
    return(none)
  }
  list(ddf = out$ddf, vcov_cov = out$vcov_theta, vcov_spcov = out$vcov_spcov, vcov_randcov = out$vcov_randcov)
}

#' Build an \code{emmeans}-compatible \code{dffun}/\code{dfargs} pair
#'
#' Used by \code{emm_basis.splm()}/\code{.spautor()} to give \code{emmeans}
#' access to Satterthwaite denominator df for arbitrary linear combinations
#' of the fixed effects -- not just the per-coefficient df already stored in
#' \code{object$ddf}, which only covers the identity contrasts. \code{emmeans}
#' calls \code{dffun(k, dfargs)} once per linear combination \code{k} (a
#' numeric vector the length of the fixed effects) it needs a denominator df
#' for, including once per row of a joint test's (QR-orthogonalized)
#' contrast matrix -- see \code{emmeans:::test.emmGrid}'s \code{joint =
#' TRUE} branch, which combines the per-row results via \code{min()}. This
#' mirrors how \code{emmeans} itself already supports Satterthwaite df for
#' \code{lme4}/\code{lmerTest} models (\code{emmeans:::emm_basis.merMod()}),
#' and is what powers \code{emmeans::joint_tests()} as well as ordinary
#' \code{emmeans()}/\code{contrast()} t-based inference.
#'
#' Mirrors \code{object$ddf}'s behavior: only attempted when
#' \code{object$ddf} is non-\code{NULL} (i.e. Satterthwaite df were
#' successfully computed when the model was fit, per its own \code{ddf}
#' argument); falls back to the package's original \code{dffun} (always
#' \code{Inf}, i.e. asymptotic/z-based inference) otherwise, or if the
#' one-time setup below fails for any reason.
#'
#' \code{emmeans::ref_grid()} always runs \code{environment(dffun) <-
#' baseenv()} on the returned \code{dffun} (to avoid retaining a
#' model-sized closure environment in the resulting \code{emmGrid} object),
#' which breaks ordinary lexical lookup of any spmodel/base function called
#' from inside \code{dffun}. To survive that, every function \code{dffun}
#' needs (\code{stats::vcov()}, \code{get_grad_g()},
#' \code{get_satterthwaite_df()}) is captured as a value in \code{dfargs}
#' up front rather than referenced by name in the closure body.
#'
#' Only \code{context} (needed for \code{get_grad_g()}'s per-contrast gradients)
#' still has to be rebuilt, since it is not itself stored on the fitted
#' object.
#'
#' @param object A fitted \code{splm}/\code{spautor} object
#'
#' @return A list with elements \code{dffun}, \code{dfargs}, and \code{mesg}
#'   (the label \code{emmeans} shows for the df method used), for direct use
#'   in the list returned by \code{emm_basis.splm()}/\code{.spautor()}
#'
#' @noRd
get_emmeans_dffun <- function(object) {
  asymptotic <- list(dffun = function(k, dfargs) Inf, dfargs = list(), mesg = "asymptotic")

  if (is.null(object$ddf)) {
    return(asymptotic)
  }

  out <- tryCatch({
    method <- get_satterthwaite_method(object, NULL)
    context <- get_cov_gradients_context(object)
    vcov_theta <- object$vcov$cov
    dfargs <- list(
      object = object, context = context, method = method, vcov_theta = vcov_theta,
      vcov_fn = stats::vcov, get_grad_g_fn = get_grad_g, get_satterthwaite_df_fn = get_satterthwaite_df
    )
    dffun <- function(k, dfargs) {
      g <- as.numeric(crossprod(k, dfargs$vcov_fn(dfargs$object)) %*% k)
      grad_g <- dfargs$get_grad_g_fn(k, dfargs$method, dfargs$context, dfargs$object)
      dfargs$get_satterthwaite_df_fn(g, grad_g, dfargs$vcov_theta)
    }
    list(dffun = dffun, dfargs = dfargs, mesg = "satterthwaite")
  }, error = function(e) NULL)

  if (is.null(out)) asymptotic else out
}
