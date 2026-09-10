#' Extract the number of observations for a fitted \code{spglm()} model object
#'
#' @description Copies \code{object$n} call from \code{nobs.splm()}
#' 
#' @param object A fitted model object from [spglm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The number of observations used for modeling. (i.e., \code{object$n}).
#'
#' @method nobs spglm
#' @noRd
#' @export
#'
#' @examples
#' binmod <- spglm(presence ~ elev, family = "binomial", data = moose, spcov_type = "exponential")
#' nobs(binmod)
nobs.spglm <- nobs.splm

#' Extract the number of observations for a fitted \code{spgautor()} model object
#'
#' @description Copies \code{object$n} call from \code{nobs.splm()}
#'
#' @param object A fitted model object from [spgautor()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The number of observations used for modeling. (i.e., \code{object$n}).
#'
#' @method nobs spgautor
#' @noRd
#' @export
#'
#' @examples
#' spgmod <- spgautor(trend ~ 1, family = "Gamma", data = seal, spcov_type = "car")
#' nobs(spgmod)
nobs.spgautor <- nobs.spautor
