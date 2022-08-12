load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_local <- FALSE # FALSE for CRAN

##### CRAN test
test_that("anova works geostatistical", {
  spmod1 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  # regular anova
  expect_error(anova(spmod1), NA)

})

#### local tests
if (test_local) {
  test_that("anova works geostatistical", {
    spmod1 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "ml")
    spmod2 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "ml")
    # regular anova
    expect_error(anova(spmod1), NA)
    expect_error(anova(spmod2), NA)
    expect_error(anova(spmod1, spmod2), NA)
    # LRT
    expect_error(tidy(anova(spmod1)), NA)
    expect_error(tidy(anova(spmod2)), NA)
    expect_error(tidy(anova(spmod1, spmod2)), NA)
    # Terms
    expect_error(anova(spmod1, Terms = c("(Intercept)", "x")), NA)
    expect_error(anova(spmod1, Terms = c(1, 2)), NA)
    all.equal(anova(spmod1, Terms = c("(Intercept)", "x")), anova(spmod1, Terms = c(1, 2)))
    # L
    L <- diag(2)
    expect_error(anova(spmod1, L = L), NA)
    all.equal(anova(spmod1, L = L)$X2, anova(spmod1, Terms = c("(Intercept)", "x"))$X2)
  })

  test_that("anova works autoregressive", {
    spmod1 <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "ml")
    spmod2 <- spautor(y ~ 1, exdata_poly, spcov_type = "car", estmethod = "ml")
    expect_error(anova(spmod1), NA)
    expect_error(anova(spmod2), NA)
    expect_error(anova(spmod1, spmod2), NA)
    expect_error(tidy(anova(spmod1)), NA)
    expect_error(tidy(anova(spmod2)), NA)
    expect_error(tidy(anova(spmod1, spmod2)), NA)
  })

  test_that("Errors appropriately return", {
    spmod1 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
    spmod2 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
    spmod3 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "ml")
    spmod4 <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls")

    # reml different fixed effects
    expect_error(anova(spmod1, spmod2))
    # reml vs ml
    expect_error(anova(spmod2, spmod3))
    # not reml and not ml
    expect_error(anova(spmod2, spmod4))
  })
}
