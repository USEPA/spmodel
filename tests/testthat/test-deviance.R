load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("Deviance works geostatistical", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_vector(deviance(spmod))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "ml")
  expect_vector(deviance(spmod))
})

test_that("Deviance works autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car")
  expect_vector(deviance(spmod))

  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "ml")
  expect_vector(deviance(spmod))
})

test_that("Deviance errors", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls")
  expect_error(deviance(spmod))
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl")
  expect_error(deviance(spmod))
})
