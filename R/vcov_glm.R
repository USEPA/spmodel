#' @rdname vcov.spmodel
#' @method vcov spglm
#' @order 3
#' @export
vcov.spglm <- vcov.splm

#' @rdname vcov.spmodel
#' @method vcov spgautor
#' @order 4
#' @export
vcov.spgautor <- vcov.spautor
