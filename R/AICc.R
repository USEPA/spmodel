#' Compute AICc of fitted model objects
#'
#' @description Compute AICc for one or
#' several fitted model objects for which a log-likelihood
#' value can be obtained.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()]
#'   where \code{estmethod} is \code{"ml"} or \code{"reml"}.
#' @param ... Optionally more fitted model objects.
#' @param k The penalty parameter, taken to be 2. Currently not allowed to differ
#'   from 2 (needed for generic consistency).
#'
#' @details When comparing models fit by maximum or restricted maximum
#'   likelihood, the smaller the AICc, the better the fit. The AICc contains
#'   a correction to AIC for small sample sizes. The theory of
#'   AICc requires that the log-likelihood has been maximized, and hence,
#'   no AICc methods exist for models where \code{estmethod} is not
#'   \code{"ml"} or \code{"reml"}. Additionally, AICc comparisons between \code{"ml"}
#'   and \code{"reml"} models are meaningless -- comparisons should only be made
#'   within a set of models estimated using \code{"ml"} or a set of models estimated
#'   using \code{"reml"}. AICc comparisons for \code{"reml"} must
#'   use the same fixed effects. To vary the covariance parameters and
#'   fixed effects simultaneously, use \code{"ml"}.
#'
#'   Hoeting et al. (2006) study AIC and AICc in a spatial context, using the AIC definition
#'   \eqn{-2loglik + 2(estparams)} and the AICc definition as
#'   \eqn{-2loglik + 2n(estparams) / (n - estparams - 1)}, where \eqn{n} is the sample size
#'   and \eqn{estparams} is the number of estimated parameters. For \code{"ml"}, \eqn{estparams} is
#'   the number of estimated covariance parameters plus the number of estimated
#'   fixed effects. For \code{"reml"}, \eqn{estparams} is the number of estimated covariance
#'   parameters.
#'
#' @return If just one object is provided, a numeric value with the corresponding
#'   AICc.
#'
#'   If multiple objects are provided, a \code{data.frame} with rows corresponding
#'   to the objects and columns representing the number of parameters estimated
#'   (\code{df}) and the AICc.
#'
#' @name AICc
#' @order 1
#' @export
#'
#' @seealso [stats::AIC()] [stats::BIC()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' AICc(spmod)
#' AIC(spmod)
#' BIC(spmod)
AICc <- function(object, ..., k = 2) {
 # method dispatch
 UseMethod("AICc", object)
}

#' @method AICc splm
#' @order 2
#' @export
AICc.splm <- function(object, ..., k = 2) {
 # set k as 2
 k <- 2
 # store object and ...
 object_list <- list(object, ...)
 # see if ... has any elements
 if (length(object_list) == 1) {
   # number of estimated parameters
   if (object$estmethod == "ml") {
     n_est_param <- object$npar + object$p
   } else {
     n_est_param <- object$npar
   }
   # error if not ml or reml
   if (!object$estmethod %in% c("ml", "reml")) {
     stop("AICc is only defined if estmethod is \"ml\" or \"reml\".", call. = FALSE)
   }
   # compute AICc
   AICc_val <- -2 * as.numeric(logLik(object)) + 2 * object$n * (n_est_param) / (object$n - n_est_param - 1)
 } else {
   # warning if ml and reml in same call
   est_methods <- vapply(object_list, function(x) x$estmethod, character(1))
   if ("ml" %in% est_methods && "reml" %in% est_methods) {
     warning("AICc and AICcc should not compare models fit with
             \"ml\" to models fit with \"reml\"", call. = FALSE)
   }
   # warning if reml and fixed effects change
   est_methods_reml <- which(est_methods == "reml")
   if (length(est_methods_reml) > 1) {
     if (any(vapply(
       est_methods_reml,
       function(x) !identical(formula(object_list[[x]]), formula(object_list[[1]])), logical(1)
     ))) {
       warning("AIC and AICc should not be used to compare models fit with \"reml\" whose fixed effect formulas differ.", call. = FALSE)
     }
   }
   # find model names provided
   object_list_names <- as.character(c(substitute(object), (as.list(substitute(list(...)))[-1])))
   # error if any names duplicated
   if (any(duplicated(object_list_names))) {
     stop("Each model object must have a unique name", call. = FALSE)
   }
   # iterate through each model
   object_AICc <- lapply(object_list, function(x) {
     # warning if estmethod not ml or reml
     if (!object$estmethod %in% c("ml", "reml")) {
       stop("AICc is only defined is estmethod is \"ml\" or \"reml\".", call. = FALSE)
     }
     if (x$estmethod == "ml") {
       n_est_param <- x$npar + x$p
     } else {
       n_est_param <- x$npar
     }
     # store degrees of freedom (parames estimated) and AICc
     data.frame(df = n_est_param, AICc = -2 * logLik(x) + 2 * x$n * (n_est_param) / (x$n - n_est_param - 1))
   })
   # put all AICc data frames together
   AICc_val <- do.call("rbind", object_AICc)
   # set rownames as model names
   row.names(AICc_val) <- object_list_names
 }
 # return AICc value
 AICc_val
}

#' @method AICc spautor
#' @order 3
#' @export
AICc.spautor <- AICc.splm
