get_delta_se <- function(fit, se.fit, family, newdata_size = 1) {
  # fit is on the link scale
  if (family %in% c("poisson", "nbinomial", "Gamma", "inverse.gaussian")) {
    g <- exp(fit)
  } else if (family %in% c("binomial", "beta")) {
    fit <- expit(fit)
    g <- fit * (1 - fit)
    # g <- exp(fit) / (1 + exp(fit))^2
  }
  val <- se.fit * g # the square root of the delta method variance

  if (family == "binomial") {
    val <- val * newdata_size
  }

  val
}
