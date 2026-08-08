#' Calculate variance-covariance matrix for a fitted model object
#'
#' @description Calculate variance-covariance matrix for a fitted model object.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param type For \code{type = "fixed"} (the default), the variance-covariance matrix
#'   of the fixed effects. If Satterthwaite degrees of freedom were calculated,
#'   \code{type = "cov"} returns the variance-covariance matrix of the covariance
#'   parameters, \code{type = "spcov"} returns the variance-covariance matrix of
#'   just the spatial variance-covariance parameters, and \code{type = "randcov"}
#'   returns the variance-covariance matrix of just the random effects (if relevant).
#' @param ... Other arguments. Not used (needed for generic consistency).
#' @param var_correct A logical indicating whether to return the corrected variance-covariance
#'   matrix for models fit using [spglm()] or [spgautor()]. The default is
#'   \code{TRUE}.
#'
#' @return The variance-covariance matrix of fixed effect 
#'   coefficients obtained via \code{coef(..., type = "fixed")}
#'   or the variance-covariance matrix of  estimated covariance parameters (when available).
#'
#' @name vcov.spmodel
#' @method vcov splm
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' vcov(spmod)
vcov.splm <- function(object, type = "fixed", ...) {
  # "cov"/"spcov"/"randcov" are only ever non-NULL when ddf = "satterthwaite"
  # succeeded at fit time -- see get_fit_ddf()/satterthwaite_core(), which
  # compute all three (the full covariance-parameter covariance matrix, and
  # its spatial-only/random-effect-only blocks) as a byproduct of the
  # denominator df themselves stored in object$ddf
  if (type == "fixed") {
    return(object$vcov$fixed)
  } else if (type == "cov") {
    return(object$vcov$cov)
  } else if (type == "spcov") {
    return(object$vcov$spcov)
  } else if (type == "randcov") {
    return(object$vcov$randcov)
  } else {
    stop("Invalid type argument. The type argument must be \"fixed\", \"cov\", \"spcov\", or \"randcov\".", call. = FALSE)
  }
}

#' @rdname vcov.spmodel
#' @method vcov spautor
#' @order 2
#' @export
vcov.spautor <- vcov.splm
