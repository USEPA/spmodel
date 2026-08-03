#' @rdname model.frame.spmodel
#' @method model.frame spglm
#' @order 3
#' @export
# spglm objects store data the same way as splm objects ($obdata), so the
# same extraction logic applies unchanged
model.frame.spglm <- model.frame.splm

#' @rdname model.frame.spmodel
#' @method model.frame spgautor
#' @order 4
#' @export
# spgautor objects store data the same way as spautor objects ($data plus
# $observed_index), so the same extraction logic applies unchanged
model.frame.spgautor <- model.frame.spautor
