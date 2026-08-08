#' @rdname cooks.distance.spmodel
#' @method cooks.distance spglm
#' @order 3
#' @export
cooks.distance.spglm <- cooks.distance.splm # GLM objects cache cooks_distance the same way as splm objects, so reuse that method

#' @rdname cooks.distance.spmodel
#' @method cooks.distance spgautor
#' @order 4
#' @export
cooks.distance.spgautor <- cooks.distance.spautor # areal GLM objects cache cooks_distance the same way as spautor objects, so reuse that method
