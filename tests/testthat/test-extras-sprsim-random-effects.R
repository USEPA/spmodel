skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

test_that("should retain intercept, slope and combined covariance", {
  data <- data.frame(cx = seq(0, 1, length.out = 12),
    cy = rep(c(0, 0.3, 0.1), 4), x = seq(-1, 1, length.out = 12),
    group = factor(rep(1:3, each = 4)), part = factor(rep(1:2, 6)))
  spatial <- spcov_params("exponential", de = 0.4, ie = 0.2, range = 0.3)
  distance <- as.matrix(dist(data[c("cx", "cy")]))
  S <- 0.4 * exp(-distance / 0.3) + diag(0.2, 12)
  same <- outer(data$group, data$group, "==")
  configurations <- list(FALSE,
    list(method_base = "base", size_base = 6, reorder_base = "none", size_new = 6, reorder_new = "none"),
    list(approximation = "vecchia", method = "all", ordering = "none"))
  components <- list(c(group = 0.6), c("x | group" = 0.6),
    c("1 | group" = 0.3, "x | group" = 0.6))
  set.seed(284)
  for (local in configurations) {
    for (j in seq_along(components)) {
      for (partition in list(NULL, ~ part)) {
        target <- S + same * switch(j, 0.6, 0.6 * tcrossprod(data$x),
          0.3 + 0.6 * tcrossprod(data$x))
        if (!is.null(partition)) target <- target * outer(data$part, data$part, "==")
        draws <- sprnorm(spatial, data = data, xcoord = cx, ycoord = cy,
          randcov_params = components[[j]], partition_factor = partition,
          local = local, samples = 6000)
        expect_equal(dim(draws), c(12L, 6000L))
        expect_lt(max(abs(rowMeans(draws))), 0.06)
        expect_lt(max(abs(cov(t(draws)) - target)), 0.09)
      }
    }
    set.seed(42)
    omitted <- sprnorm(spatial, data = data, xcoord = cx, ycoord = cy, local = local)
    set.seed(42)
    explicit <- sprnorm(spatial, data = data, xcoord = cx, ycoord = cy, local = local, randcov_params = NULL)
    expect_identical(explicit, omitted)
  }
  set.seed(91)
  latent <- sprnorm(spatial, data = data, xcoord = cx, ycoord = cy, randcov_params = c(group = 0.6))
  expected <- rpois(12, exp(latent))
  set.seed(91)
  actual <- sprpois(spatial, data = data, xcoord = cx, ycoord = cy, randcov_params = c(group = 0.6))
  expect_equal(as.numeric(actual), expected)
})

test_that("direct covariance accepts NULL and slope-only components", {
  data <- data.frame(x = seq(-1, 1, length.out = 12), group = factor(rep(1:3, 4)))
  for (spatial in list(spcov_params("none", ie = 0.2), spcov_params("car", de = 0.4, ie = 0.2, range = 0.2, extra = 0))) {
    W <- matrix(1, 12, 12) - diag(12)
    set.seed(7)
    omitted <- sprnorm(spatial, data = data, W = W)
    set.seed(7)
    explicit <- sprnorm(spatial, data = data, W = W, randcov_params = NULL)
    expect_identical(omitted, explicit)
    set.seed(8)
    baseline <- sprnorm(spatial, data = data, W = W, samples = 6000)
    set.seed(9)
    slope <- sprnorm(spatial, data = data, W = W, samples = 6000, randcov_params = c("x | group" = 0.6))
    target <- 0.6 * tcrossprod(data$x) * outer(data$group, data$group, "==")
    expect_lt(max(abs(cov(t(slope)) - cov(t(baseline)) - target)), 0.09)
  }
})

test_that("generator formulas retain their local environment and reject ambiguous components", {
  local_env <- environment()
  formula <- get_simulation_random(c("1 | group" = 0.3, "x | group" = 0.6), env = local_env)
  expect_identical(environment(formula), local_env)
  expect_equal(get_randcov_names(formula), c("1 | group", "x | group"))
  expect_error(get_simulation_random(0.3, local_env), "name for each")
  expect_error(get_simulation_random(c(group = 0.3, group = 0.6), local_env), "No /")
  expect_error(get_simulation_random(c("group/subgroup" = 0.3), local_env), "only one variance")
})
