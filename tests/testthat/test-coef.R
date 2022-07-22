load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_local <- FALSE # FALSE for CRAN

##### CRAN test
test_that("coef works geostatistical", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_error(coef(spmod), NA)
})

#### local tests
if (test_local) {
test_that("coef works geostatistical", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_error(coef(spmod), NA)
    expect_error(coef(spmod, type = "spcov"), NA)
    expect_null(coef(spmod, type = "randcov"))
    expect_error(coefficients(spmod), NA)
    expect_error(coefficients(spmod, type = "spcov"), NA)
    expect_null(coefficients(spmod, type = "randcov"))


    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, random = ~group)
    expect_error(coef(spmod), NA)
    expect_error(coef(spmod, type = "spcov"), NA)
    expect_error(coef(spmod, type = "randcov"), NA)
    expect_error(coefficients(spmod), NA)
    expect_error(coefficients(spmod, type = "spcov"), NA)
    expect_error(coefficients(spmod, type = "randcov"), NA)
  })

  test_that("coef works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_error(coef(spmod), NA)
    expect_error(coef(spmod, type = "spcov"), NA)
    expect_null(coef(spmod, type = "randcov"))
    expect_error(coefficients(spmod), NA)
    expect_error(coefficients(spmod, type = "spcov"), NA)
    expect_null(coefficients(spmod, type = "randcov"))

    spmod <- spautor(y ~ x, exdata_poly, "car", random = ~group)
    expect_error(coef(spmod), NA)
    expect_error(coef(spmod, type = "spcov"), NA)
    expect_error(coef(spmod, type = "randcov"), NA)
    expect_error(coefficients(spmod), NA)
    expect_error(coefficients(spmod, type = "spcov"), NA)
    expect_error(coefficients(spmod, type = "randcov"), NA)
  })

  test_that("errors return", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_error(coef(spmod, type = "xyz"))
    expect_error(coefficients(spmod, type = "xyz"))
  })
}
