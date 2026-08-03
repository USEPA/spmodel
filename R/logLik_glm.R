#' @rdname logLik.spmodel
#' @method logLik spglm
#' @order 3
#' @export
# the log-likelihood is extracted from object$optim$value the same way
# regardless of family, so the GLM methods reuse the Gaussian implementations
logLik.spglm <- logLik.splm

#' @rdname logLik.spmodel
#' @method logLik spgautor
#' @order 4
#' @export
logLik.spgautor <- logLik.spautor
