#' @rdname confint.spmodel
#' @method confint spglm
#' @order 3
#' @export
# spglm/spgautor reuse the splm/spautor Wald interval logic; intervals are
# on the link scale since coef()/vcov() return link-scale fixed effects
confint.spglm <- confint.splm

#' @rdname confint.spmodel
#' @method confint spgautor
#' @order 4
#' @export
confint.spgautor <- confint.spautor
