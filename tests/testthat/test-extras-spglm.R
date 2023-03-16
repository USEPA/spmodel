test_that("blank test", {
  expect_null(NULL)
})

test_local <- FALSE # FALSE for CRAN

if (test_local) {

  set.seed(1)

  # SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  # add variables
  n <- NROW(exdata)
  exdata$bern <- rbinom(n, size = 1, prob = 0.5)
  exdata$size <- 5
  exdata$bin <- rbinom(n, size = exdata$size, prob = 0.5)
  exdata$prop <- runif(n)
  exdata$count <- rpois(n, lambda = 5)
  exdata$cont <- rgamma(n, shape = 1, rate = 1)

  #### CRAN checks
  test_that("the model runs for binomial data", {
    expect_error(spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
    expect_error(spglm(cbind(bin, size) ~ x, family = "binomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "ml"), NA)
  })

  test_that("the model runs for proportion data", {
    expect_error(spglm(prop ~ x, family = "beta", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
    expect_error(spglm(prop ~ x, family = beta, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"), NA)
  })

  test_that("the model runs for count data", {
    expect_error(spglm(count ~ x, family = poisson, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml"), NA)
    expect_error(spglm(count ~ x, family = "nbinomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"), NA)
  })

  test_that("the model runs for continuous data", {
    expect_error(spglm(cont ~ x, family = "Gamma", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml"), NA)
    expect_error(spglm(cont ~ x, family = inverse.gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"), NA)
    expect_error(spglm(cont ~ x, family = gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"), NA)
  })

  test_that("the model runs on other data sets", {
    expect_error(spglm(abs(y) ~ x, family = "Gamma", data = exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
    expect_error(spglm(abs(y) ~ x, family = Gamma, data = exdata_poly, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml"), NA)

  })



}
