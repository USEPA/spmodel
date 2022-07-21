load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("loocv works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_vector(loocv(spmod))
  expect_vector(loocv(spmod, local = TRUE))
  # cores 2 for cran check
  expect_vector(loocv(spmod, local = list(parallel = TRUE, ncores = 2)))
  expect_equal(length(loocv(spmod, cv_predict = TRUE)), 2)
  expect_equal(length(loocv(spmod, cv_predict = TRUE, local = TRUE)), 2)
  expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE)), 3)
  expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = TRUE)), 3)
  # cores 2 for cran check
  expect_equal(length(loocv(spmod, cv_predict = TRUE, local = list(parallel = TRUE, ncores = 2))), 2)
  expect_equal(length(loocv(spmod, cv_predict = TRUE, local = list(parallel = TRUE, method = "all", ncores = 2))), 2)
  expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = list(parallel = TRUE, ncores = 2))), 3)
  expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = list(parallel = TRUE, method = "all", ncores = 2))), 3)

  # random effects
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, random = ~group)
  expect_vector(loocv(spmod))
  expect_vector(loocv(spmod, local = TRUE))
})

test_that("loocv works auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_vector(loocv(spmod))
  # cores 2 for cran check
  expect_vector(loocv(spmod, local = list(parallel = TRUE, ncores = 2)))
  expect_equal(length(loocv(spmod, cv_predict = TRUE)), 2)
  expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE)), 3)
  # cores 2 for cran check
  expect_equal(length(loocv(spmod, cv_predict = TRUE, local = list(parallel = TRUE, ncores = 2))), 2)
  expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = list(parallel = TRUE, ncores = 2))), 3)

  # random effects
  spmod <- spautor(y ~ x, exdata_poly, "car", random = ~group)
  expect_vector(loocv(spmod))
  expect_vector(loocv(spmod, local = TRUE))
})
