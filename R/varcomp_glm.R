# spglm/spgautor reuse the splm/spautor varcomp logic directly (same variance
# decomposition applies on the link-function scale)
#' @rdname varcomp
#' @method varcomp spglm
#' @order 4
#' @export
varcomp.spglm <- varcomp.splm

#' @rdname varcomp
#' @method varcomp spgautor
#' @order 5
#' @export
varcomp.spgautor <- varcomp.spautor
