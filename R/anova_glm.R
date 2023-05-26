#' @rdname anova.spmodel
#' @method anova spglm
#' @order 3
#' @export
anova.spglm <- anova.splm

#' @rdname anova.spmodel
#' @method anova spgautor
#' @order 4
#' @export
anova.spgautor <- anova.spautor

#' @rdname anova.spmodel
#' @method tidy anova.spglm
#' @order 7
#' @export
tidy.anova.spglm <- tidy.anova.splm

#' @rdname anova.spmodel
#' @method tidy anova.spgautor
#' @order 8
#' @export
tidy.anova.spgautor <- tidy.anova.spautor
