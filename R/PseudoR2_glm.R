#' @rdname pseudoR2
#' @method pseudoR2 spglm
#' @order 4
#' @export
# spglm/spgautor reuse the splm/spautor pseudo r-squared formula since both
# store a precomputed object$pseudoR2 (deviance ratio) the same way
pseudoR2.spglm <- pseudoR2.splm

#' @rdname pseudoR2
#' @method pseudoR2 spgautor
#' @order 5
#' @export
pseudoR2.spgautor <- pseudoR2.spautor
