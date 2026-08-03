#' @rdname anova.spmodel
#' @method anova spglm
#' @order 3
#' @export
# spglm/spgautor reuse the splm/spautor anova() and tidy() logic directly,
# since the GLHT and LRT computations only rely on generics (vcov, logLik,
# coefficients, model.matrix) that are defined for all model classes
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
