load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("AIC and AICc works geostatistical", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_vector(AIC(spmod))
  expect_vector(AICc(spmod))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "ml")
  expect_vector(AIC(spmod))
  expect_vector(AICc(spmod))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls")
  expect_error(AIC(spmod))
  expect_error(AICc(spmod))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl")
  expect_error(AIC(spmod))
  expect_error(AICc(spmod))
})

test_that("AIC and AICc works autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")
  expect_vector(AIC(spmod))
  expect_vector(AICc(spmod))

  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "ml")
  expect_vector(AIC(spmod))
  expect_vector(AICc(spmod))
})


test_that("AIC and AICc works for two models", {
  spmod0 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  spmod1 <- splm(y ~ 1, exdata, spcov_type = "none", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_error(AIC(spmod0, spmod1), NA)
  expect_error(AICc(spmod0, spmod1), NA)
})


test_that("Errors appropriately return", {
  spmod1 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")

  spmod0 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls")
  expect_error(AIC(spmod0))
  expect_error(AICc(spmod0))
  expect_error(AIC(spmod0, spmod1))
  expect_error(AICc(spmod0, spmod1))

  spmod0 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl")
  expect_error(AIC(spmod0))
  expect_error(AICc(spmod0))
  expect_error(AIC(spmod0, spmod1))
  expect_error(AICc(spmod0, spmod1))

  spmod0 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  spmod1 <- splm(y ~ 1, exdata, spcov_type = "none", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_error(AIC(spmod0, spmod1, spmod0))
  expect_error(AICc(spmod0, spmod1, spmod0))
})

test_that("Warnings appropriately return", {
  spmod0 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  spmod1 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  spmod2 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "ml")
  expect_warning(AIC(spmod0, spmod1))
  expect_warning(AIC(spmod1, spmod2))
  expect_warning(AICc(spmod0, spmod1))
  expect_warning(AICc(spmod1, spmod2))
})
