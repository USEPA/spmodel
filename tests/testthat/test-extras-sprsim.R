skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

test_that("sprnorm() local$approximation = 'vecchia' works", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  spcov_params_val <- spcov_params("exponential", de = 1, ie = 0.05, range = 1.5)

  sim1 <- sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "vecchia", size = 10))
  expect_true(inherits(sim1, "matrix"))
  expect_equal(dim(sim1), c(NROW(exdata), 20))
  expect_true(all(is.finite(sim1)))

  # method = "all" (no truncation) should closely match the exact (local =
  # FALSE) unconditional distribution -- see get_sprnorm_vecchia()
  R <- 4000
  set.seed(1)
  sim_exact <- sprnorm(spcov_params_val, samples = R, data = exdata, xcoord = xcoord, ycoord = ycoord, local = FALSE)
  set.seed(2)
  sim_vecchia_all <- sprnorm(spcov_params_val, samples = R, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "vecchia", method = "all"))
  expect_equal(rowMeans(sim_exact), rowMeans(sim_vecchia_all), tolerance = 0.1)
  expect_equal(apply(sim_exact, 1, sd), apply(sim_vecchia_all, 1, sd), tolerance = 0.15)

  # neighbor-selection rules and distance/covariance truncation both run
  expect_vector(sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "vecchia", size = 5, method = "distance"))[, 1])
  expect_vector(sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "vecchia", size = 5, method = "covariance"))[, 1])

  # random effects supported via covmatrix() reuse
  exdata_re <- exdata
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  randcov_params_val <- randcov_params(grp = 0.3)
  sim_re <- sprnorm(spcov_params_val, samples = 20, data = exdata_re, xcoord = xcoord, ycoord = ycoord, randcov_params = randcov_params_val, local = list(approximation = "vecchia", size = 10))
  expect_true(all(is.finite(sim_re)))

  # size_base/size_new/reorder_base/kmeans_new (renamed from reorder/kmeans)
  # all still work under the default ("low-rank") approximation
  expect_vector(sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(size_base = 50, size_new = 20, reorder_base = "random", kmeans_new = FALSE))[, 1])

  # invalid local$approximation errors informatively
  expect_error(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "bogus")), "local\\$approximation must be")
})