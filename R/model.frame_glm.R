#' @rdname model.frame.spmodel
#' @method model.frame spglm
#' @order 3
#' @export
model.frame.spglm <- model.frame.splm

#' @rdname model.frame.spmodel
#' @method model.frame spgautor
#' @order 4
#' @export
model.frame.spgautor <- model.frame.spautor
