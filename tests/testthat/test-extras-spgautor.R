test_that("blank test", {
  expect_null(NULL)
})

test_local <- FALSE # FALSE for CRAN

if (test_local) {

  set.seed(1)

  # SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

  # add variables
  n <- NROW(exdata_poly)
  exdata_poly$bern <- rbinom(n, size = 1, prob = 0.5)
  exdata_poly$size <- 5
  exdata_poly$bin <- rbinom(n, size = exdata_poly$size, prob = 0.5)
  exdata_poly$prop <- runif(n)
  exdata_poly$count <- rpois(n, lambda = 5)
  exdata_poly$cont <- rgamma(n, shape = 1, rate = 1)

  #### CRAN checks
  test_that("the model runs for binomial data", {
    expect_error(spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml"), NA)
    expect_error(spgautor(cbind(bin, size) ~ x, family = "binomial", data = exdata_poly, spcov_type = "sar", estmethod = "ml"), NA)
  })

  test_that("the model runs for proportion data", {
    expect_error(spgautor(prop ~ x, family = "beta", data = exdata_poly, spcov_type = "car", estmethod = "reml"), NA)
    expect_error(spgautor(prop ~ x, family = beta, data = exdata_poly, spcov_type = "sar", estmethod = "reml"), NA)
  })

  test_that("the model runs for count data", {
    expect_error(spgautor(count ~ x, family = poisson, data = exdata_poly, spcov_type = "sar", estmethod = "reml"), NA)
    expect_error(spgautor(count ~ x, family = "nbinomial", data = exdata_poly, spcov_type = "car", estmethod = "ml"), NA)
  })

  test_that("the model runs for continuous data", {
    expect_error(spgautor(cont ~ x, family = "Gamma", data = exdata_poly, spcov_type = "car", estmethod = "reml"), NA)
    expect_error(spgautor(cont ~ x, family = inverse.gaussian, data = exdata_poly, spcov_type = "sar", estmethod = "ml"), NA)
    expect_error(spgautor(cont ~ x, family = gaussian, data = exdata_poly, spcov_type = "car", estmethod = "reml"), NA)
  })

  test_that("the model runs on other data sets", {
    expect_error(spgautor(abs(y) ~ x, family = "Gamma", data = exdata_Mpoly, spcov_type = "car", estmethod = "reml"), NA)
    expect_error(spgautor(abs(y) ~ x, family = Gamma, data = exdata_Upoly, spcov_type = "sar", estmethod = "reml"), NA)

  })



}
