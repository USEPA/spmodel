#' @method AICc spglm
#' @order 4
#' @export
# spglm and spgautor reuse the splm/spautor AICc computation directly since
# npar and p are stored the same way across all four model classes
AICc.spglm <- AICc.splm

#' @method AICc spgautor
#' @order 5
#' @export
AICc.spgautor <- AICc.spautor
