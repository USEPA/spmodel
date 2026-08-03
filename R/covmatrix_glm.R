#' @rdname covmatrix
#' @method covmatrix spglm
#' @order 3
#' @export
covmatrix.spglm <- covmatrix.splm # GLM objects have the same structure needed here (obdata, spcov params, etc.), so the splm method is reused as-is

#' @rdname covmatrix
#' @method covmatrix spgautor
#' @order 3
#' @export
covmatrix.spgautor <- covmatrix.spautor # areal GLM objects have the same structure needed here, so the spautor method is reused as-is
