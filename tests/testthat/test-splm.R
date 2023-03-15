set.seed(1)

# SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

#### CRAN checks
test_that("the model runs for several scenarios", {
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "gaussian", estmethod = "sv-cl"), NA)
})
