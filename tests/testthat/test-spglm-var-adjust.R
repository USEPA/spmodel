g8_spglm_fixture <- function() {
  set.seed(1808)
  data <- data.frame(cx = runif(96), cy = runif(96), x = rnorm(96))
  eta <- 0.7 + 0.35 * data$x + sprnorm(spcov_params("exponential", de = .7, ie = .1, range = .3),
    data = data, xcoord = cx, ycoord = cy, local = FALSE)
  data$count <- rpois(nrow(data), exp(eta))
  data$binary <- rbinom(nrow(data), 1, plogis(eta - .7))
  list(data = data, index = 1L + (data$cx >= .5) + 2L * (data$cy >= .5))
}

test_that("local spglm propagates latent variance independently of var_adjust", {
  fixture <- g8_spglm_fixture()
  methods <- c("none", "theoretical", "pooled", "empirical")
  for (family in c("poisson", "binomial")) {
    formula <- if (family == "poisson") count ~ x else binary ~ x
    fits <- lapply(methods, function(method) {
      spglm(formula, family = as.character(family), data = fixture$data, xcoord = cx, ycoord = cy,
        spcov_type = "exponential", local = list(index = fixture$index, var_adjust = method))
    })
    ref <- fits[[1L]]
    latent <- vcov(ref) - vcov(ref, var_correct = FALSE)
    for (fit in fits[-1L]) {
      expect_identical(coef(fit), coef(ref))
      expect_identical(coef(fit, "spcov"), coef(ref, "spcov"))
      expect_equal(vcov(fit) - vcov(fit, var_correct = FALSE), latent, tolerance = 1e-10)
    }
  }
})
