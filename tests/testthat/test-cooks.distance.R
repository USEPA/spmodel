load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("Cook's distance works geostatistical", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_vector(cooks.distance(spmod))
  expect_true(min(cooks.distance(spmod)) >= 0)
})

test_that("Cook's distance works autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car")
  expect_vector(cooks.distance(spmod))
  expect_true(min(cooks.distance(spmod)) >= 0)
})
