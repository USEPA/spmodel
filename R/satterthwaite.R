#' Compute Satterthwaite denominator degrees of freedom
#'
#' @description Compute Satterthwaite denominator degrees of freedom
#'   \eqn{t}-based (rather than asymptotic \eqn{z}-based)
#'   fixed effect inference in small samples.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Satterthwaite degrees of freedom are generally more appropriate than
#'   asymptotic degrees of freedom for small samples. They can be computationally costly
#'   for sample sizes exceeding 500; however, for sample sizes this large, they Satterthwaite
#'   and asymptotic degrees of freedom should yield very similar inferences.
#'
#' @return A named numeric vector of Satterthwaite degrees of freedom for each
#'   fixed effect.
#'
#' @seealso [splm()] [spautor()] [anova.spmodel()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y, estmethod = "reml"
#' )
#' satterthwaite(spmod)
#' }
#' 
#' @references
#'   Rencher, Alvin C. and Schaalje, G. Bruce (2008). Linear Models in 
#'   Statistics, Second Edition. John Wiley & Sons.
satterthwaite <- function(object, ...) {
  UseMethod("satterthwaite")
}

#' @rdname satterthwaite
#' @param method The method by which to compute gradients. \code{"numeric"}
#'   for numerical differentiation and \code{"closed"} for closed form solutions.
#'   The default \code{"closed"} for \code{"exponential"}, \code{"gaussian"},
#'   \code{"spherical"}, \code{"none"}, and \code{"ie"} spatial covariance
#'   functions (without anisotropy) and \code{"numeric"} otherwise.
#' @method satterthwaite splm
#' @order 2
#' @export
satterthwaite.splm <- function(object, method, ...) {
  if (missing(method)) method <- NULL
  satterthwaite_core(object, method)$ddf
}

#' Shared Satterthwaite computation for \code{satterthwaite()} and \code{get_fit_ddf()}
#'
#' Does the actual work behind \code{satterthwaite.splm()}/\code{.spautor()},
#' factored out so \code{get_fit_ddf()} (called by \code{splm()}/
#' \code{spautor()} at fit time) can also get at \code{vcov_theta} -- which it
#' stores in \code{object$vcov$cov} (see [vcov.spmodel()]) -- without paying
#' for the expensive context/Hessian computation a second time.
#'
#' @param object A fitted \code{splm}/\code{spautor} object
#' @param method See \code{satterthwaite()}
#'
#' @return A list with elements \code{ddf} (as returned by
#'   \code{satterthwaite()}) and \code{vcov_theta} (the estimated covariance
#'   matrix of the free covariance parameters, or \code{NULL} if it could not
#'   be computed)
#'
#' @noRd
satterthwaite_core <- function(object, method) {
  method <- get_satterthwaite_method(object, method)

  validate_satterthwaite_scope(object, method)
  # context bundles every quantity needed to treat the fitted covariance
  # matrix Sigma as a function of its free parameters theta (spatial + random
  # effect variance components) -- built once and reused below since it does
  # not depend on which coefficient/contrast is being tested
  context <- get_cov_gradients_context(object)

  # Cov(theta_hat): the delta-method machinery below needs the sampling
  # uncertainty of the estimated covariance parameters themselves, not just
  # of the fixed effects -- asymptotic intervals completetly ignore this uncertainty
  # while Satterthwaite's df estimate shrink when this uncertainty is large
  vcov_theta <- get_vcov_theta(method, context, object)

  betahat <- coef(object, type = "fixed")
  p <- length(betahat)

  # one row per fixed effect coefficient (an identity contrast Li picks out
  # coefficient i alone); satterthwaite_anova()/fai_cornelius() instead build
  # more general (possibly multi-row) L matrices for joint hypotheses
  L <- diag(p)
  colnames(L) <- names(betahat)
  rownames(L) <- names(betahat)

  ddf <- lapply(seq_len(p), function(i) {
    Li <- L[i, ]
    # g = Var(Li'betahat) = the squared standard error already reported by
    # summary()/vcov() -- Satterthwaite treats g_hat as (approximately)
    # g * chisq_df / df for some unknown df, and estimates df by matching the
    # variance of that chi-squared surrogate to the delta-method variance of
    # g_hat computed just below
    g <- as.numeric(crossprod(Li, vcov(object)) %*% Li)
    # gradient of g(theta) with respect to the free covariance parameters --
    # the delta method turns uncertainty in thetahat (vcov_theta) into
    # uncertainty in g_hat via Var(g_hat) ~= grad_g' Cov(theta_hat) grad_g
    grad_g <- get_grad_g(Li, method, context, object)
    satterthwaite_df <- get_satterthwaite_df(g, grad_g, vcov_theta)
  })

  ddf <- unlist(ddf)
  names(ddf) <- names(betahat)

  # split the joint covariance-parameter covariance matrix into its spatial
  # (spcov) and random effect (randcov) blocks, for vcov(object, type =
  # "spcov"/"randcov") -- context$cov_names_free is exactly the spcov names
  # followed by the randcov names (see get_cov_gradients_context()), so
  # simple name-based subsetting recovers each block; NULL propagates
  # through whenever vcov_theta itself could not be computed
  if (is.null(vcov_theta)) {
    vcov_spcov <- NULL
    vcov_randcov <- NULL
  } else {
    vcov_spcov <- if (length(context$spcov_names_free)) {
      vcov_theta[context$spcov_names_free, context$spcov_names_free, drop = FALSE]
    } else {
      NULL
    }
    vcov_randcov <- if (length(context$randcov_names_free)) {
      vcov_theta[context$randcov_names_free, context$randcov_names_free, drop = FALSE]
    } else {
      NULL
    }
  }

  list(ddf = ddf, vcov_theta = vcov_theta, vcov_spcov = vcov_spcov, vcov_randcov = vcov_randcov)
}

#' @rdname satterthwaite
#' @method satterthwaite spautor
#' @order 3
#' @export
satterthwaite.spautor <- satterthwaite.splm

# "closed" and "numeric" only differ in how Cov(theta_hat) and grad_g get
# computed (get_vcov_theta()/get_grad_g()); the Satterthwaite formula itself
# (get_satterthwaite_df()) is identical either way. Closed-form derivatives
# of Sigma are only implemented for the covariance types listed below without
# anisotropy (the dSig_dtheta_spcov.<type>() methods in get_dSig_dtheta.R);
# every other covariance type/anisotropic fit falls back to numerically
# differentiating the likelihood/covariance itself.

get_satterthwaite_method <- function(object, method) {

  satterthwaite_closed_form_types <- c("exponential", "gaussian", "spherical", "none", "ie")

  # anisotropy makes the distance matrix itself a function of theta (rotate/
  # scale), which the closed-form dSig_dtheta_spcov.<type>() derivatives do
  # not account for -- so anisotropic fits must fall back to "numeric" even
  # when the covariance type itself has a closed form
  has_closed_form <- inherits(coef(object, type = "spcov"), satterthwaite_closed_form_types) &&
    !isTRUE(object$anisotropy)
  if (is.null(method)) {
    if (has_closed_form) method <-"closed" else method <-"numeric" 
  }
  if (method == "closed" && !has_closed_form) {
    method <- "numeric"
    warning("Closed form not available for this model. Using method = \"numeric\".", call. = FALSE)
  }
  method
}

# Satterthwaite's approximation: treat g_hat (a variance estimate) as
# approximately g * chisq_df / df for some unknown df. Under that
# approximation Var(g_hat) = 2 * g^2 / df, so matching Var(g_hat) to its
# delta-method estimate (grad_g' Cov(theta_hat) grad_g) and solving for df
# gives the formula below. A larger delta-method variance (more uncertain
# covariance parameters relative to g) implies a smaller estimated df, i.e.
# a heavier-tailed t reference distribution.
get_satterthwaite_df <- function(g, grad_g, vcov_theta) {
  # vcov_theta is NULL whenever Cov(theta_hat) could not be computed (a
  # non-positive-definite Hessian/Fisher information, already warned about
  # where vcov_theta is built) -- crossprod(grad_g, NULL) does not error, it
  # silently returns a length-mismatched result, so this must be checked
  # explicitly rather than left to propagate
  if (is.null(vcov_theta)) {
    return(NA)
  }
  as.numeric(2 * g^2 / (crossprod(grad_g, vcov_theta) %*% grad_g))
}