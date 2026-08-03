#' Delta-method standard error for a response-scale prediction
#'
#' @param fit A fitted value on the link scale
#' @param se.fit The link-scale standard error
#' @param family The response family
#' @param newdata_size The binomial trial size (used only when \code{family} is \code{"binomial"})
#'
#' @return The (approximate) response-scale standard error, obtained by scaling
#'   the link-scale standard error by the derivative of the inverse link function
#'
#' @noRd
get_delta_se <- function(fit, se.fit, family, newdata_size = 1) {
  # fit is on the link scale (w)
  # delta method: Var(g(fit)) ~ g'(fit)^2 * Var(fit), so the response-scale SE
  # is se.fit * |g'(fit)|; g here is the derivative of the inverse link
  # function for each family (log link -> exp(fit);
  # logit/logistic-type -> fit(1-fit) = p(1 - p) = exp(fit)/(1 + exp(fit))^2)
  if (family %in% c("poisson", "nbinomial", "Gamma", "inverse.gaussian")) {
    g <- exp(fit)
  } else if (family %in% c("binomial", "beta")) {
    fit <- expit(fit)
    g <- fit * (1 - fit)
  }
  val <- se.fit * g # the square root of the delta method variance

  # binomial response is a proportion on the link scale; rescale back up to
  # the count scale by multiplying by the number of trials
  if (family == "binomial") {
    val <- val * newdata_size
  }

  val
}
