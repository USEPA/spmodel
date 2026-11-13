skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

test_that("simulation skips unused coefficient covariance adjustments", {
  local_mocked_bindings(get_W_ij = function(...) stop("unused coefficient adjustment"))
  set.seed(66)
  d <- data.frame(cx = runif(240), cy = runif(240))
  params <- spcov_params("exponential", de = .7, ie = .3, range = .3)
  {
    expect_length(sprnorm(params, data = d, xcoord = cx, ycoord = cy), 240)
  }
})

test_that("iid simulation separates slopes and partition-specific effects without dense algebra", {
  d <- data.frame(x = seq(-1, 1, length.out = 24), g = factor(rep(1:3, 8)), part = factor(rep(1:2, each = 12)))
  local_mocked_bindings(cov_matrix = function(...) stop("dense covariance"),
    partition_matrix = function(...) stop("dense partition"),
    splm = function(...) stop("unnecessary model fit"))
  target <- (.4 + .3 * tcrossprod(d$x)) * outer(d$g, d$g, "==") * outer(d$part, d$part, "==") + diag(.2, 24)
  parameters <- list(spcov_params("none", ie = .2), spcov_params("ie", ie = .2),
    spcov_params("exponential", de = 0, ie = .2, range = .3))
  reference <- NULL
  for (params in parameters) {
    set.seed(714)
    draws <- sprnorm(params, data = d, randcov_params = c("1 | g" = .4, "x | g" = .3),
      partition_factor = ~ part, samples = 15000)
    if (is.null(reference)) reference <- draws else expect_identical(draws, reference)
    expect_lt(max(abs(cov(t(draws)) - target)), .05)
  }
  large <- data.frame(g = factor(seq_len(20000)), part = factor(rep(1:2, 10000)))
  expect_length(sprnorm(parameters[[1]], data = large, randcov_params = c(g = .4),
    partition_factor = ~ part), 20000)
  zero <- spcov_params("none", ie = 0)
  expect_identical(sprnorm(zero, data = d, mean = 2), rep(2, 24))
  expect_true(all(is.finite(sprnorm(zero, data = d, randcov_params = c(g = .4)))))
  bad <- d
  bad$g <- as.numeric(bad$g)
  expect_error(sprnorm(parameters[[1]], data = bad, randcov_params = c(g = .4)), "categorical or factor")
})
