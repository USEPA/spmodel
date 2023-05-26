#' @rdname model.matrix.spmodel
#' @method model.matrix spglm
#' @order 3
#' @export
model.matrix.spglm <- model.matrix.splm

#' @rdname model.matrix.spmodel
#' @method model.matrix spgautor
#' @order 4
#' @export
model.matrix.spgautor <- model.matrix.spautor
