#' @rdname model.matrix.spmodel
#' @method model.matrix spglm
#' @order 3
#' @export
# same underlying data storage as splm, so the same construction logic applies
model.matrix.spglm <- model.matrix.splm

#' @rdname model.matrix.spmodel
#' @method model.matrix spgautor
#' @order 4
#' @export
# same underlying data storage as spautor, so the same construction logic applies
model.matrix.spgautor <- model.matrix.spautor
