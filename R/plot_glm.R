#' @rdname plot.spmodel
#' @method plot spglm
#' @order 3
#' @export
# diagnostic plots only use generics (fitted(), rstandard(), etc.), which
# already dispatch correctly on spglm objects, so the splm plotting code
# works unchanged
plot.spglm <- plot.splm

#' @rdname plot.spmodel
#' @method plot spgautor
#' @order 4
#' @export
# same reasoning as plot.spglm above, but for the spautor/spgautor pairing
plot.spgautor <- plot.spautor
