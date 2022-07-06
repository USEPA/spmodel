load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("confint works geostatistical", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_error(confint(spmod), NA)
  expect_error(confint(spmod, parm = "(Intercept)", level = 0.99), NA)
  expect_error(confint(spmod, parm = "x", level = 0.90), NA)
})

test_that("confint works autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_error(confint(spmod), NA)
  expect_error(confint(spmod, parm = "(Intercept)", level = 0.99), NA)
  expect_error(confint(spmod, parm = "x", level = 0.90), NA)
})
