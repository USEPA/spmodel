#' @rdname deviance.spmodel
#' @method deviance spglm
#' @order 3
#' @export
deviance.spglm <- deviance.splm # glm variant reuses the same reml/ml check and lookup as the linear case

#' @rdname deviance.spmodel
#' @method deviance spgautor
#' @order 4
#' @export
deviance.spgautor <- deviance.spautor
