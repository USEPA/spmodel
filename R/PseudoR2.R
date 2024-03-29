#' Compute a pseudo r-squared
#'
#' @description Compute a pseudo r-squared for a fitted model object.
#'
#' @param object A fitted model object from [splm()], [spautor()], [spglm()], or [spgautor()].
#' @param adjust A logical indicating whether the pseudo r-squared
#'   should be adjusted to account for the number of explanatory variables. The
#'   default is \code{FALSE}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Several pseudo r-squared statistics exist for in the literature.
#'   We define this pseudo r-squared as one minus the ratio of the deviance of a full model
#'   relative to the deviance of a null (intercept only) model. This pseudo r-squared
#'   can be viewed as a generalization of the classical r-squared definition
#'   seen as one minus the ratio of error sums of squares from the full model relative
#'   to the error sums of squares from the null model. If adjusted, the adjustment
#'   is analogous to the the classical r-squared adjustment.
#'
#' @return The pseudo r-squared as a numeric vector.
#'
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' pseudoR2(spmod)
pseudoR2 <- function(object, ...) {
  UseMethod("pseudoR2", object)
}

#' @rdname pseudoR2
#' @method pseudoR2 splm
#' @order 2
#' @export
pseudoR2.splm <- function(object, adjust = FALSE, ...) {
  if (adjust) {
    has_intercept <- "(Intercept)" %in% tidy(object)$term
    pr2 <- object$pseudoR2
    pr2_adj <- 1 - (1 - pr2) * (object$n - 1 * has_intercept) / (object$n - object$p)
    return(pr2_adj)
  } else {
    return(object$pseudoR2)
  }
}

#' @rdname pseudoR2
#' @method pseudoR2 spautor
#' @order 3
#' @export
pseudoR2.spautor <- pseudoR2.splm
