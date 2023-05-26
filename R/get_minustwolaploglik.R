get_minustwolaploglik <- function(laploglik_products, estmethod, n, p, spcov_profiled = FALSE, randcov_profiled = FALSE) {
  # spcov profiled and randcov profiled always FALSE right now
  if (estmethod == "reml") {
    minustwolaploglik <- as.numeric(laploglik_products$l00 + laploglik_products$l01 + laploglik_products$l1 +
      laploglik_products$l2 + laploglik_products$l3 +
      (n - p) * log(2 * pi))
  } else if (estmethod == "ml") {
    minustwolaploglik <- as.numeric(laploglik_products$l00 + laploglik_products$l01 + laploglik_products$l1 +
      laploglik_products$l2 +
      n * log(2 * pi))
  }

  minustwolaploglik
}
