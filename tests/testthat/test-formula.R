load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("formula works geostatistical", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_equal(formula(spmod), y ~ x)

  spmod <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_equal(formula(spmod), y ~ 1)
})

test_that("formula works autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car")
  expect_equal(formula(spmod), y ~ x)

  spmod <- spautor(y ~ 1, exdata_poly, spcov_type = "car")
  expect_equal(formula(spmod), y ~ 1)
})
