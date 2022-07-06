load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("logLik works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_vector(logLik(spmod))
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, estmethod = "ml")
  expect_vector(logLik(spmod))
})

test_that("logLik works auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_vector(logLik(spmod))
  spmod <- spautor(y ~ x, exdata_poly, "car", estmethod = "ml")
  expect_vector(logLik(spmod))
})

test_that("logLik appropriately errors", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, estmethod = "sv-wls")
  expect_error(logLik(spmod))
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, estmethod = "sv-cl")
  expect_error(logLik(spmod))
})
