#' @rdname formula.spmodel
#' @method formula spglm
#' @order 3
#' @export
# GLM-type models store formula identically to splm, so reuse that method
formula.spglm <- formula.splm

#' @rdname formula.spmodel
#' @method formula spgautor
#' @order 4
#' @export
# reuse the spautor method since spgautor objects store formula the same way
formula.spgautor <- formula.spautor
