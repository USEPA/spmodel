#' @rdname covmatrix
#' @method covmatrix spglm
#' @order 3
#' @export
covmatrix.spglm <- covmatrix.splm

#' @rdname covmatrix
#' @method covmatrix spgautor
#' @order 3
#' @export
covmatrix.spgautor <- covmatrix.spautor
