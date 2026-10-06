test_that("splm conditional simulation returns requested outputs", {
  data <- data.frame(cx = seq_len(12), x = rep(c(-1, 0, 1), 4),
    y = c(2, 1, 4, 3, 6, 2, 1, 4, 5, 2, 3, 6))
  fit <- splm(y ~ x, data, xcoord = cx, ddf = "asymptotic",
    spcov_initial = spcov_initial("exponential", de = 0.4, ie = 0.2, range = 3, known = "given"))
  draws <- conditional(fit, data[1:2, ], samples = 2, output = "all")
  expect_named(draws, c("newdata", "beta", "object"))
  expect_equal(dim(draws$newdata), c(2, 2))
  expect_true(all(is.finite(draws$newdata)))
})

test_that("Poisson conditional simulation supports exact and local paths", {
  data <- data.frame(cx = seq_len(12), x = rep(c(-1, 0, 1), 4),
    y = c(2, 1, 4, 3, 6, 2, 1, 4, 5, 2, 3, 6))
  fit <- spglm(y ~ x, data, family = "poisson", xcoord = cx,
    spcov_initial = spcov_initial("exponential", de = 0.4, ie = 0.2, range = 3, known = "given"))
  for (local in list(FALSE, list(approximation = "low-rank", size_base = 6, reorder_base = "none"),
    list(approximation = "vecchia", size = 3, ordering = "none"))) {
    draws <- conditional(fit, data[1:2, ], samples = 2, type = "new", local = local)
    expect_equal(dim(draws), c(2, 2))
    expect_true(all(is.finite(draws) & draws >= 0 & draws == floor(draws)))
  }
})

test_that("Areal conditional simulation uses fitted missing-response locations", {
  data <- data.frame(x = seq(-1, 1, length.out = 12),
    y = c(2, NA, 4, 3, 6, 2, 1, 4, NA, 2, 3, 6))
  W <- 1 * (abs(outer(seq_len(12), seq_len(12), "-")) == 1)
  covariance <- spcov_initial("car", de = 0.4, ie = 0.2, range = 0.2, known = "given")
  fits <- list(spautor(y ~ x, data, W = W, spcov_initial = covariance, ddf = "asymptotic"),
    spgautor(y ~ x, data, W = W, family = "poisson", spcov_initial = covariance))
  for (fit in fits) {
    draws <- conditional(fit, samples = 2)
    expect_equal(dim(draws), c(2, 2))
    expect_identical(rownames(draws), as.character(fit$missing_index))
    expect_true(all(is.finite(draws)))
  }
})
