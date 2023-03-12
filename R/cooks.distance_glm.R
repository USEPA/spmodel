#' @rdname cooks.distance.spmodel
#' @method cooks.distance spglm
#' @order 3
#' @export
cooks.distance.spglm <- cooks.distance.splm

#' @rdname cooks.distance.spmodel
#' @method cooks.distance spgautor
#' @order 4
#' @export
cooks.distance.spgautor <- cooks.distance.spautor
