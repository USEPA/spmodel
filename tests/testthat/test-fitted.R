load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_local <- FALSE # FALSE for CRAN

##### CRAN test
test_that("fitted distance works geostatistical", {

  # no random effect
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_vector(fitted(spmod))
})

#### local tests
if (test_local) {
test_that("fitted distance works geostatistical", {

  # no random effect
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_vector(fitted(spmod))
  expect_vector(fitted.values(spmod))
  expect_vector(fitted(spmod, type = "spcov"))
  expect_vector(fitted.values(spmod, type = "spcov"))
  expect_null(fitted(spmod, type = "randcov"))
  expect_null(fitted.values(spmod, type = "randcov"))

  # random effect
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~group)
  expect_vector(fitted(spmod))
  expect_vector(fitted.values(spmod))
  expect_vector(fitted(spmod, type = "spcov"))
  expect_vector(fitted.values(spmod, type = "spcov"))
  expect_vector(fitted(spmod, type = "randcov"))
  expect_vector(fitted.values(spmod, type = "randcov"))
})

test_that("fitted distance works autoregressive", {

  # no random effect
  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car")
  expect_vector(fitted(spmod))
  expect_vector(fitted.values(spmod))
  expect_vector(fitted(spmod, type = "spcov"))
  expect_vector(fitted.values(spmod, type = "spcov"))
  expect_null(fitted(spmod, type = "randcov"))
  expect_null(fitted.values(spmod, type = "randcov"))

  # random effect
  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", random = ~group)
  expect_vector(fitted(spmod))
  expect_vector(fitted.values(spmod))
  expect_vector(fitted(spmod, type = "spcov"))
  expect_vector(fitted.values(spmod, type = "spcov"))
  expect_vector(fitted(spmod, type = "randcov"))
  expect_vector(fitted.values(spmod, type = "randcov"))
})

test_that("errors return", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_error(fitted(spmod, type = "xyz"))
  expect_error(fitted.values(spmod, type = "xyz"))
})
}
