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
  exdata$bernfac <- factor(ifelse(exdata$bern == 0, "a", "b"))
  exdata$size <- 5
  exdata$bin <- rbinom(n, size = exdata$size, prob = 0.5)
  exdata$prop <- runif(n)
  exdata$count <- rpois(n, lambda = 5)
  exdata$cont <- rgamma(n, shape = 1, rate = 1)
  exdata$offset <- 1.2

  test_that("the model runs for binomial data", {
    spgmod <- spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
    expect_s3_class(spgmod, "spglm")
    expect_vector(AUC(spgmod))
    spgmod <- spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
    expect_s3_class(spgmod, "spglm")
    expect_vector(AUC(spgmod))
    expect_error(spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
    spgmod <- spglm(cbind(bin, size) ~ x, family = "binomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "ml")
    expect_s3_class(spgmod, "spglm")
    expect_error(AUC(spgmod))
    spgmod <- spglm(y > 0 ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
    expect_s3_class(spgmod, "spglm")
    expect_vector(AUC(spgmod))


    # complicated models
    expect_error(spglm(bern ~ x + offset(offset), family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE), NA)
    # check when data have an unequal number of local observations
    expect_error(spglm(bern ~ x + offset(offset), family = binomial, data = exdata[-1, , drop = FALSE], xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE), NA)
    expect_error(spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE), NA)
    expect_error(spglm(cbind(bin, size) ~ x, family = "binomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "ml", partition_factor = ~group), NA)
  })

  test_that("the model runs for proportion data", {
    expect_error(spglm(prop ~ x, family = "beta", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
    expect_error(spglm(prop ~ x, family = beta, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"), NA)

    # complicated models
    expect_error(spglm(prop ~ x, family = "beta", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", anisotropy = TRUE), NA)
    expect_error(spglm(prop ~ x + offset(offset),
      family = beta, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml",
      local = list(method = "kmeans")
    ), NA)
  })

  test_that("the model runs for count data", {
    spgmod <- spglm(count ~ x, family = poisson, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml")
    expect_s3_class(spgmod, "spglm")
    expect_error(AUC(spgmod))
    expect_error(spglm(count ~ x, family = "nbinomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"), NA)

    # complicated models
    expect_error(spglm(count ~ x + offset(offset),
      family = poisson, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml",
      partition_factor = ~group, local = TRUE
    ), NA)
    expect_error(spglm(count ~ x, family = "nbinomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml", random = ~group), NA)
  })

  test_that("the model runs for continuous data", {
    expect_error(spglm(cont ~ x, family = "Gamma", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml"), NA)
    expect_error(spglm(cont ~ x, family = inverse.gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"), NA)

    # SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
    expect_error(spglm(cont ~ x, family = gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"))


    # complicated models
    expect_error(spglm(cont ~ x + offset(offset),
      family = "Gamma", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml",
      random = ~ group + subgroup, local = TRUE
    ), NA)
    expect_error(spglm(cont ~ x, family = inverse.gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml", partition_factor = ~subgroup, anisotropy = TRUE), NA)
    ## SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
    expect_error(spglm(cont ~ x,
      family = gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml",
      random = ~group
    ))
  })

  test_that("the model runs on other data sets", {
    expect_error(spglm(abs(y) ~ x, family = "Gamma", data = exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
    expect_error(spglm(abs(y) ~ x, family = Gamma, data = exdata_poly, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml"), NA)


    # complicated models
    expect_error(spglm(abs(y) ~ x, family = "Gamma", data = exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, partition_factor = ~group, anisotropy = TRUE, local = TRUE), NA)
    expect_error(spglm(abs(y) ~ x,
      family = Gamma, data = exdata_poly, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml",
      random = ~group, partition_factor = ~subgroup
    ), NA)
  })
}
