#' @rdname labels.spmodel
#' @method labels spglm
#' @order 3
#' @export
# fixed-effect term labels come from the formula alone, so the GLM methods
# just reuse the Gaussian model implementations
labels.spglm <- labels.splm

#' @rdname labels.spmodel
#' @method labels spgautor
#' @order 4
#' @export
labels.spgautor <- labels.spautor
