#' Glance at a fitted model object
#'
#' @description Returns a row of model
#'   summaries from a fitted model object. Glance returns the same number of columns for all models
#'   and estimation methods. If a particular summary is undefined for a model
#'   or estimation method (e.g., likelihood statistics for estimation methods
#'   \code{"sv-wls"} or \code{"sv-cl"}) of \code{splm()} objects), \code{NA} is returned for that summary.
#'
#' @param x A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A single-row tibble with columns
#'   \itemize{
#'     \item{\code{n}}{ The sample size.}
#'     \item{\code{p}}{ The number of fixed effects.}
#'     \item{\code{npar}}{ The number of estimated covariance parameters.}
#'     \item{\code{value}}{ The optimized value of the fitting function}
#'     \item{\code{AIC}}{ The AIC.}
#'     \item{\code{AICc}}{ The AICc.}
#'     \item{\code{logLik}}{ The log-likelihood}
#'     \item{\code{deviance}}{ The deviance.}
#'     \item{\code{pseudo.r.squared}}{ The pseudo r-squared}
#'   }
#'
#' @name glance.spmodel
#' @method glance splm
#' @order 1
#' @export
#'
#' @seealso [AIC.spmodel()] [AICc()] [logLik.spmodel()] [deviance.spmodel()] [pseudoR2()] [tidy.spmodel()] [augment.spmodel()]
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' glance(spmod)
glance.splm <- function(x, ...) {
  is_likbased <- x$estmethod %in% c("ml", "reml")
  tibble::tibble(
    n = x$n,
    p = x$p,
    npar = x$npar,
    value = x$optim$value,
    AIC = ifelse(is_likbased, AIC(x), NA),
    AICc = ifelse(is_likbased, AICc(x), NA),
    logLik = ifelse(is_likbased, logLik(x), NA),
    deviance = ifelse(is_likbased, deviance(x), NA),
    pseudo.r.squared = pseudoR2(x),
    # cv.crit = loocv(x)
  )
}

#' @rdname glance.spmodel
#' @method glance spautor
#' @order 2
#' @export
glance.spautor <- glance.splm
