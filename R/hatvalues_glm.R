#' @rdname hatvalues.spmodel
#' @method hatvalues spglm
#' @order 3
#' @export
# GLM hat values are stored the same way as the Gaussian models, so the
# splm/spautor extractors are reused as-is
hatvalues.spglm <- hatvalues.splm

#' @rdname hatvalues.spmodel
#' @method hatvalues spgautor
#' @order 4
#' @export
hatvalues.spgautor <- hatvalues.spautor
