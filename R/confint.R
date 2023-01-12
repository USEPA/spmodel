#' Confidence intervals for fitted model parameters
#'
#' @description Computes confidence intervals for one or more parameters in a fitted
#'   model object.
#'
#' @param object A fitted model object from [splm()] or [spautor()].
#' @param parm A specification of which parameters are to be given confidence
#'   intervals (a character vector of names). If missing, all parameters are considered.
#' @param level The confidence level required. The default is \code{0.95}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return Gaussian-based confidence intervals (two-sided and equal-tailed) for the
#'   fixed effect coefficients based on the confidence level specified by \code{level}.
#'
#' @method confint splm
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' confint(spmod)
#' confint(spmod, parm = "waterY", level = 0.90)
confint.splm <- function(object, parm, level = 0.95, ...) {
  # if (type == "fixed") ## may add spcov and randcov confidence intervals later
  alpha <- 1 - level
  # tstar <- qt(1 - alpha / 2, df = object$n - object$p)
  tstar <- qnorm(1 - alpha / 2)
  estimates <- coef(object, type = "fixed")
  variances <- diag(vcov(object, type = "fixed"))
  lower <- estimates - tstar * sqrt(variances)
  upper <- estimates + tstar * sqrt(variances)
  confints <- cbind(lower, upper)
  rownames(confints) <- names(estimates)
  colnames(confints) <- c(paste(alpha / 2 * 100, "%"), paste((1 - alpha / 2) * 100, "%"))
  if (missing(parm)) {
    return(confints)
  } else {
    return(confints[row.names(confints) %in% parm, , drop = FALSE])
  }
}

#' @rdname confint.splm
#' @method confint spautor
#' @export
confint.spautor <- confint.splm
