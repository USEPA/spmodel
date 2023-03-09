#' @rdname anova.spmodel
#' @method anova spglm
#' @export
anova.spglm <- anova.splm

#' @rdname anova.spmodel
#' @method anova spgautor
#' @export
anova.spgautor <- anova.spautor

#' @rdname anova.spmodel
#' @method tidy anova.spglm
#' @export
tidy.anova.spglm <- tidy.anova.splm

#' @rdname anova.spmodel
#' @method tidy anova.spgautor
#' @export
tidy.anova.spgautor <- tidy.anova.spautor
