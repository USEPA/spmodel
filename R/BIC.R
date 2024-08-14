# #' Compute BIC of fitted model objects
# #'
# #' @description Compute BIC for one or
# #' several fitted model objects for which a log-likelihood
# #' value can be obtained.
# #'
# #' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()]
# #'   where \code{estmethod} is \code{"ml"} or \code{"reml"}.
# #' @param ... Optionally more fitted model objects.
# #'
# #' @details When comparing models fit by maximum or restricted maximum
# #'   likelihood, the smaller the BIC, the better the fit. The theory of
# #'   BIC requires that the log-likelihood has been maximized, and hence,
# #'   no BIC methods exist for models where \code{estmethod} is not
# #'   \code{"ml"} or \code{"reml"}. Additionally, BIC comparisons between \code{"ml"}
# #'   and \code{"reml"} models are meaningless -- comparisons should only be made
# #'   within a set of models estimated using \code{"ml"} or a set of models estimated
# #'   using \code{"reml"}. BIC comparisons for \code{"reml"} must
# #'   use the same fixed effects. To vary the covariance parameters and
# #'   fixed effects simultaneously, use \code{"ml"}.
# #'
# #'   BIC is defined as \eqn{-2loglik + log(n)(estparams)}, where \eqn{n} is the sample size
# #'   and \eqn{estparams} is the number of estimated parameters. For \code{"ml"}, \eqn{estparams} is
# #'   the number of estimated covariance parameters plus the number of estimated
# #'   fixed effects. For \code{"reml"}, \eqn{estparams} is the number of estimated covariance
# #'   parameters.
# #'
# #' @return If just one object is provided, a numeric value with the corresponding
# #'   BIC.
# #'
# #'   If multiple objects are provided, a \code{data.frame} with rows corresponding
# #'   to the objects and columns representing the number of parameters estimated
# #'   (\code{df}) and the BIC.
# #'
# #' @name BIC.spmodel
# #' @method BIC splm
# #' @order 1
# #' @export
# #'
# #' @examples
# #' spmod <- splm(z ~ water + tarp,
# #'   data = caribou,
# #'   spcov_type = "exponential", xcoord = x, ycoord = y
# #' )
# #' BIC(spmod)
# BIC.splm <- function(object, ...) {
#
#   # set k as log of sample size
#   k <- log(object$n)
#
#   # store object and ...
#   object_list <- list(object, ...)
#
#   # see if ... has any elements
#   if (length(object_list) == 1) {
#
#     # number of estimated parameters
#     if (object$estmethod == "ml") {
#       n_est_param <- object$npar + object$p
#     } else {
#       n_est_param <- object$npar
#     }
#
#     # error if not ml or reml
#     if (!object$estmethod %in% c("ml", "reml")) {
#       stop("BIC is only defined is estmethod is \"ml\" or \"reml\".", call. = FALSE)
#     }
#     # compute BIC
#     BIC_val <- -2 * logLik(object) + k * (n_est_param)
#   } else {
#
#
#     # warning if ml and reml in same call
#     est_methods <- vapply(object_list, function(x) x$estmethod, character(1))
#     if ("ml" %in% est_methods && "reml" %in% est_methods) {
#       warning("BIC and BICc should not compare models fit with
#               \"ml\" to models fit with \"reml\"", call. = FALSE)
#     }
#
#     # warning if reml and fixed effects change
#     est_methods_reml <- which(est_methods == "reml")
#     if (length(est_methods_reml) > 1) {
#       if (any(vapply(
#         est_methods_reml,
#         function(x) !identical(formula(object_list[[x]]), formula(object_list[[1]])), logical(1)
#       ))) {
#         warning("BIC and BICc should not be used to compare models fit with \"reml\" whose fixed effect formulas differ.", call. = FALSE)
#       }
#     }
#
#     # find model names provided
#     object_list_names <- as.character(c(substitute(object), (as.list(substitute(list(...)))[-1])))
#
#     # error if any names duplicated
#     if (any(duplicated(object_list_names))) {
#       stop("Each model object must have a unique name", call. = FALSE)
#     }
#     # iterate through each model
#     object_BIC <- lapply(object_list, function(x) {
#       # warning if estmethod not ml or reml
#       if (!object$estmethod %in% c("ml", "reml")) {
#         stop("BIC is only defined is estmethod is \"ml\" or \"reml\".", call. = FALSE)
#       }
#
#       if (x$estmethod == "ml") {
#         n_est_param <- x$npar + x$p
#       } else {
#         n_est_param <- x$npar
#       }
#
#       # store degrees of freedom (parames estimated) and BIC
#       data.frame(df = n_est_param, BIC = -2 * logLik(x) + k * (n_est_param))
#     })
#     # put all BIC data frames together
#     BIC_val <- do.call("rbind", object_BIC)
#     # set rownames as model names
#     row.names(BIC_val) <- object_list_names
#   }
#   # return BIC value
#   BIC_val
# }
#
# #' @rdname BIC.spmodel
# #' @method BIC spautor
# #' @order 2
# #' @export
# BIC.spautor <- BIC.splm
#
