satterthwaite_anova <- function(object, ...) {
  UseMethod("satterthwaite_anova")
}

#' @noRd
#' @exportS3Method
# generalizes satterthwaite()'s per-coefficient df to (possibly multi-row,
# i.e. multi-degree-of-freedom) joint hypotheses: a Wald chi-squared
# statistic on q df, divided by q, is already an F(q, Inf) statistic, so
# reusing get_marginal_Chi2()'s existing Chi2 and just swapping the reference
# distribution's denominator df from Inf to a Satterthwaite/Fai-Cornelius
# estimate (DenDF, from fai_cornelius() below) turns the same statistic into
# a small-sample-adjusted F test
satterthwaite_anova.splm <- function(object, ..., test = TRUE, Terms, L, method) {

  if (missing(Terms)) Terms <- NULL
  if (missing(L)) L <- NULL
  if (missing(method)) method <- NULL

  method <- get_satterthwaite_method(object, method)

  validate_satterthwaite_scope(object, method)

  # one hypothesis matrix per term (or per explicit Terms/L request) -- see
  # get_L(); each element of L is tested/df-adjusted independently below
  L <- get_L(L, Terms, object)

  context <- get_cov_gradients_context(object)
  vcov_theta <- get_vcov_theta(method, context, object)
  betahat <- coef(object, type = "fixed")

  DenDF <- unlist(lapply(L, fai_cornelius, betahat = betahat,
                 vcov_theta = vcov_theta, context = context,
                 method = method, object))

  # call anova() function and appropriate scale quantities to be F instead of X2
  anova_val <- do.call(rbind, lapply(L, get_marginal_Chi2, object))
  anova_val$`F value`<- anova_val$Chi2/anova_val$Df
  anova_val$NumDF <- anova_val$Df
  anova_val$DenDF <- DenDF
  anova_val <- anova_val[, c("NumDF", "DenDF", "F value")]
  anova_val$`Pr(>F)` <- pf(anova_val$`F value`, anova_val$NumDF, anova_val$DenDF, lower.tail = FALSE)
  anova_val
}

# get_cov_gradients_context()/get_vcov_theta()/fai_cornelius() (and, inside
# fai_cornelius(), get_grad_g()) all dispatch on object's class themselves
# so method is reused as-is rather than duplicated
#' @noRd
#' @exportS3Method
satterthwaite_anova.spautor <- satterthwaite_anova.splm

fai_cornelius <- function(Lsub, betahat, vcov_theta, context, method, object) {
  UseMethod("fai_cornelius", object)
}

#' @noRd
#' @exportS3Method
# Fai & Cornelius (1996): Satterthwaite's df formula only applies to a single
# (scalar) variance estimate, but a joint test on q > 1 rows of Lsub needs
# one overall denominator df for the whole F test. The approach first
# rotates Lsub's rows into q *uncorrelated* linear combinations eta_i (an
# eigendecomposition of C = Var(Lsub %*% betahat) diagonalizes their
# covariance), so each eta_i gets its own valid single-row Satterthwaite df
# nu_i exactly as in satterthwaite(). Those q individual df estimates are
# then combined into one DenDF via Fai & Cornelius's formula, which behaves
# like a (weighted harmonic) mean -- E is a stabilizing sum of nu_i/(nu_i-2)
# terms, and DenDF = 2E/(E-q) reduces to nu_1 itself when q = 1.
fai_cornelius.splm <- function(Lsub, betahat, vcov_theta, context, method, object) {

  if (!is.matrix(Lsub)) Lsub <- matrix(Lsub, nrow = 1)
  q <- nrow(Lsub)
  C <- Lsub %*% tcrossprod(vcov(object), Lsub)
  eig <- eigen(C)
  U <- eig$vectors
  eta <- crossprod(U, Lsub)

  nu <- lapply(seq_len(q), function(i) {
    eta_i <- eta[i, ]
    g <- as.numeric(crossprod(eta_i, vcov(object)) %*% eta_i)
    grad_g <- get_grad_g(eta_i, method, context, object)
    satterthwaite_df <- get_satterthwaite_df(g, grad_g, vcov_theta)
  })

  nu <- unlist(nu)

  # the combination formula is only valid when every component's own df
  # exceeds 2 (nu_i/(nu_i - 2), a component's contribution to E, is only
  # finite/positive there) -- a component's own Satterthwaite df collapsing
  # to <= 2 signals its variance is too poorly estimated for this joint test
  # to be valid, so DenDF (and the resulting F test) is reported
  # as NA rather than silently producing a (potentially) misleading statistic.
  if (any(!is.finite(nu)) || any(nu <= 2)) {
    DenDF <- NA
  } else {
    E <- sum(nu / (nu - 2))
    DenDF <- if (E > q) 2 * E / (E - q) else NA
  }
  DenDF
}

# get_grad_g() dispatches on object's class themselves 
# so method is reused as-is rather than duplicated
#' @noRd
#' @exportS3Method
fai_cornelius.spautor <- fai_cornelius.splm
