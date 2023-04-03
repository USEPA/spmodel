test_that("blank test", {
  expect_null(NULL)
})

test_local <- FALSE # FALSE for CRAN

if (test_local) {

  # SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

  ##############################################################################
  ############################ AIC and AICc (test-aic.R)
  ##############################################################################

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

  test_that("Matches for lm", {
    spmod <- splm(y ~ x, exdata, "none", estmethod = "ml")
    lmod <- lm(y ~ x, exdata)
    expect_equal(AIC(spmod), AIC(lmod))
  })

  ##############################################################################
  ############################ anova (test-anova.R)
  ##############################################################################

  spmod1 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  # regular anova
  expect_error(anova(spmod1), NA)

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

  ##############################################################################
  ############################ augment (test-augment.R)
  ##############################################################################

  test_that("augment works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_error(augment(spmod), NA)
    expect_error(augment(spmod, se_fit = TRUE), NA)
    expect_s3_class(augment(spmod), "tbl")
    expect_error(augment(spmod, newdata = newexdata), NA)
    augment(spmod, newdata = newexdata, interval = "confidence")
    augment(spmod, newdata = newexdata, interval = "prediction")
    expect_error(augment(spmod, newdata = newexdata, se_fit = TRUE), NA)
    augment(spmod, newdata = newexdata, interval = "confidence", se_fit = TRUE)
    augment(spmod, newdata = newexdata, interval = "prediction", se_fit = TRUE)
    expect_s3_class(augment(spmod, newdata = newexdata), "tbl")
  })

  test_that("augment works auto", {
    spmod <- spautor(y ~ x, exdata_Mpoly, "car")
    expect_error(augment(spmod), NA)
    expect_error(augment(spmod, se_fit = TRUE), NA)
    expect_s3_class(augment(spmod), "tbl")
    expect_s3_class(augment(spmod), "sf")
    expect_error(augment(spmod, newdata = spmod$newdata), NA)
    augment(spmod, newdata = spmod$newdata, interval = "confidence")
    augment(spmod, newdata = spmod$newdata, interval = "prediction")
    expect_error(augment(spmod, newdata = spmod$newdata, se_fit = TRUE), NA)
    augment(spmod, newdata = spmod$newdata, interval = "confidence", se_fit = TRUE)
    augment(spmod, newdata = spmod$newdata, interval = "prediction", se_fit = TRUE)
    expect_s3_class(augment(spmod, newdata = spmod$newdata), "tbl")
    expect_s3_class(augment(spmod, newdata = spmod$newdata), "sf")
  })

  test_that("augment works with drop", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    aug_mod <- augment(spmod) # default drop = TRUE
    expect_equal(NCOL(aug_mod), 9)
    aug_mod <- augment(spmod, drop = FALSE)
    expect_equal(NCOL(aug_mod), 11)

    aug_pred <- augment(spmod, newdata = newexdata)
    expect_equal(NCOL(aug_pred), 6) # default drop = FALSE

    spmod <- spautor(y ~ x, exdata_Mpoly, "car")
    aug_mod <- augment(spmod) # default drop = TRUE
    expect_equal(NCOL(aug_mod), 8)
    aug_mod <- augment(spmod, drop = FALSE)
    expect_equal(NCOL(aug_mod), 11)

    aug_pred <- augment(spmod, newdata = spmod$newdata)
    expect_equal(NCOL(aug_pred), 7) # default drop = FALSE
  })

  test_that("augment works with types", {
    exdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"))
    newexdata_sf <- sf::st_as_sf(newexdata, coords = c("xcoord", "ycoord"))

    # newdata output type same as data output type

    # df fit df pred
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    aug_mod <- augment(spmod) # default drop = TRUE
    expect_true(inherits(aug_mod, "tbl"))
    expect_false(inherits(aug_mod, "sf"))
    aug_mod <- augment(spmod, drop = FALSE)
    expect_true(inherits(aug_mod, "tbl"))
    expect_false(inherits(aug_mod, "sf"))
    aug_pred <- augment(spmod, newdata = newexdata)
    expect_true(inherits(aug_pred, "tbl"))
    expect_false(inherits(aug_pred, "sf"))

    # sf fit sf pred
    spmod <- splm(y ~ x, exdata_sf, "exponential", xcoord, ycoord)
    aug_mod <- augment(spmod) # default drop = TRUE
    expect_true(inherits(aug_mod, "tbl"))
    expect_true(inherits(aug_mod, "sf"))
    aug_mod <- augment(spmod, drop = FALSE)
    expect_true(inherits(aug_mod, "tbl"))
    expect_true(inherits(aug_mod, "sf"))
    aug_pred <- augment(spmod, newdata = newexdata_sf)
    expect_true(inherits(aug_pred, "tbl"))
    expect_true(inherits(aug_pred, "sf"))

    # df fit sf pred
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    aug_mod <- augment(spmod) # default drop = TRUE
    expect_true(inherits(aug_mod, "tbl"))
    expect_false(inherits(aug_mod, "sf"))
    aug_mod <- augment(spmod, drop = FALSE)
    expect_true(inherits(aug_mod, "tbl"))
    expect_false(inherits(aug_mod, "sf"))
    aug_pred <- augment(spmod, newdata = newexdata_sf)
    expect_true(inherits(aug_pred, "tbl"))
    expect_false(inherits(aug_pred, "sf"))

    # sf fit df pred
    newexdata$.xcoord <- newexdata$xcoord
    newexdata$.ycoord <- newexdata$ycoord
    spmod <- splm(y ~ x, exdata_sf, "exponential", xcoord, ycoord)
    aug_mod <- augment(spmod) # default drop = TRUE
    expect_true(inherits(aug_mod, "tbl"))
    expect_true(inherits(aug_mod, "sf"))
    aug_mod <- augment(spmod, drop = FALSE)
    expect_true(inherits(aug_mod, "tbl"))
    expect_true(inherits(aug_mod, "sf"))
    aug_pred <- augment(spmod, newdata = newexdata)
    expect_true(inherits(aug_pred, "tbl"))
    expect_true(inherits(aug_pred, "sf"))



    # sf in
    spmod <- spautor(y ~ x, exdata_Mpoly, "car")
    expect_true(inherits(aug_mod, "tbl"))
    expect_true(inherits(aug_mod, "sf"))
    aug_mod <- augment(spmod, drop = FALSE)
    expect_true(inherits(aug_mod, "tbl"))
    expect_true(inherits(aug_mod, "sf"))
    aug_pred <- augment(spmod, newdata = spmod$newexdata)
    expect_true(inherits(aug_pred, "tbl"))
    expect_true(inherits(aug_pred, "sf"))

    # df in
    W <- 1 * sf::st_intersects(exdata_Mpoly, sparse = FALSE)
    diag(W) <- 0
    exdata_Mpoly_df <- st_drop_geometry(exdata_Mpoly)
    spmod <- spautor(y ~ x, exdata_Mpoly_df, "car", W = W)
    aug_mod <- augment(spmod)
    expect_true(inherits(aug_mod, "tbl"))
    expect_false(inherits(aug_mod, "sf"))
    aug_mod <- augment(spmod, drop = FALSE)
    expect_true(inherits(aug_mod, "tbl"))
    expect_false(inherits(aug_mod, "sf"))
    aug_pred <- augment(spmod, newdata = spmod$newdata)
    expect_true(inherits(aug_pred, "tbl"))
    expect_false(inherits(aug_pred, "sf"))
  })

  ##############################################################################
  ############################ coef (test-coef.R)
  ##############################################################################

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

  ##############################################################################
  ############################ confint (test-confint.R)
  ##############################################################################

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

  ##############################################################################
  ############################ cooks distance (test-cooks.distance.R)
  ##############################################################################

  test_that("Cook's distance works geostatistical", {
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
    expect_vector(cooks.distance(spmod))
    expect_true(min(cooks.distance(spmod)) >= 0)
  })

  test_that("Cook's distance works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car")
    expect_vector(cooks.distance(spmod))
    expect_true(min(cooks.distance(spmod)) >= 0)
  })

  ##############################################################################
  ############################ covmatrix (test-covmatrix.R)
  ##############################################################################

  test_that("covmatrix works geostatistical", {
    spcov_type <- "exponential"
    spmod <- splm(y ~ x, exdata_M, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, random = ~ group)
    expect_equal(c(99, 99), dim(covmatrix(spmod)))
    expect_equal(c(1, 99), dim(covmatrix(spmod, newdata = spmod$newdata)))
    expect_equal(c(10, 99), dim(covmatrix(spmod, newdata = newexdata)))
  })

  test_that("covmatrix works autoregressive", {
    spcov_type <- "car"
    spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type, random = ~ group)
    expect_equal(c(48, 48), dim(covmatrix(spmod)))
    expect_equal(c(1, 48), dim(covmatrix(spmod, newdata = spmod$newdata)))
  })

  test_that("covmatrix works geostatistical", {
    spcov_type <- "exponential"
    spmod <- splm(y ~ x, exdata_M, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, partition_factor = ~ group)
    expect_equal(c(99, 99), dim(covmatrix(spmod)))
    expect_equal(c(1, 99), dim(covmatrix(spmod, newdata = spmod$newdata)))
    expect_equal(c(10, 99), dim(covmatrix(spmod, newdata = newexdata)))
  })

  test_that("covmatrix works autoregressive", {
    spcov_type <- "car"
    spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type, partition_factor = ~ group)
    expect_equal(c(48, 48), dim(covmatrix(spmod)))
    expect_equal(c(1, 48), dim(covmatrix(spmod, newdata = spmod$newdata)))
  })

  ##############################################################################
  ############################ deviance (test-deviance.R)
  ##############################################################################

  test_that("Deviance works geostatistical", {
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
    expect_vector(deviance(spmod))

    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "ml")
    expect_vector(deviance(spmod))
  })

  test_that("Deviance works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car")
    expect_vector(deviance(spmod))

    spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "ml")
    expect_vector(deviance(spmod))
  })

  test_that("Deviance errors", {
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls")
    expect_error(deviance(spmod))
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl")
    expect_error(deviance(spmod))
  })

  ##############################################################################
  ############################ esv (test-esv.R)
  ##############################################################################

  test_that("esv works", {

    # regular implementation
    esv1 <- esv(y ~ x, exdata, xcoord, ycoord)
    expect_s3_class(esv1, "data.frame")
    expect_equal(NROW(esv1), 15)
    expect_equal(NCOL(esv1), 4)

    esv1_q <- esv(y ~ x, exdata, "xcoord", "ycoord")
    expect_s3_class(esv1_q, "data.frame")
    expect_equal(NROW(esv1_q), 15)
    expect_equal(NCOL(esv1_q), 4)

    # quoting works
    expect_equal(esv1, esv1_q)

    # specifying bins and cutoff
    esv2 <- esv(y ~ x, exdata, xcoord, ycoord, bins = 30, cutoff = 5)
    expect_s3_class(esv2, "data.frame")
    expect_equal(NROW(esv2), 30)
    expect_equal(NCOL(esv1), 4)
    dist_matrix <- spdist(exdata, "xcoord", "ycoord")

    # specifying distance matrix
    esv3 <- esv(y ~ x, exdata, dist_matrix = dist_matrix)
    expect_s3_class(esv3, "data.frame")
    expect_equal(NROW(esv3), 15)
    expect_equal(NCOL(esv1), 4)

    # specifying partition factor
    esv4 <- esv(y ~ x, exdata, xcoord, ycoord, partition_factor = ~group)
    expect_s3_class(esv1, "data.frame")
    expect_equal(NROW(esv1), 15)
    expect_equal(NCOL(esv1), 4)
    expect_false(identical(esv1, esv4)) # make sure results are not identical to full esv

    # works with sf object
    exdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"))
    expect_error(esv(y ~ x, exdata_sf), NA)

    # works with one dimension
    expect_error(esv(y ~ x, exdata, xcoord), NA)

    # errors occur
    expect_error(esv(y ~ x, exdata))
    expect_error(esv(y ~ x, exdata, xcoord_xyz))
    expect_error(esv(y ~ x, exdata, xcoord, ycoord_xyz))
  })

  ##############################################################################
  ############################ fitted (test-fitted.R)
  ##############################################################################

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

  ##############################################################################
  ############################ formula (test-formula.R)
  ##############################################################################

  test_that("formula works geostatistical", {
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
    expect_equal(formula(spmod), y ~ x)

    spmod <- splm(y ~ 1, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
    expect_equal(formula(spmod), y ~ 1)
  })

  test_that("formula works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car")
    expect_equal(formula(spmod), y ~ x)

    spmod <- spautor(y ~ 1, exdata_poly, spcov_type = "car")
    expect_equal(formula(spmod), y ~ 1)
  })

  ##############################################################################
  ############################ FREE GENERICS (test-generics.R)
  ##############################################################################

  test_that("free generics work geostatistical", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_true(is.language(terms(spmod))) # works because there is an object spmod$terms
    expect_true(is.call(getCall(spmod))) # works because there is an object spmod$call
  })

  test_that("free generics work geostatistical", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_true(is.language(terms(spmod))) # works because there is an object spmod$terms
    expect_true(is.call(getCall(spmod))) # works because there is an object spmod$call
    expect_equal(update(formula(spmod), y ~ 1), y ~ 1) # works because formula(object) works (object has class spmod)
    spmod <- update(spmod, y ~ 1, spcov_type = "spherical")
    expect_s3_class(spmod, "splm")
    expect_equal(formula(spmod), y ~ 1)
    expect_s3_class(coefficients(spmod, type = "spcov"), "spherical")
    spmod <- update(spmod, . ~ . + offset(x))
    expect_vector(model.offset(model.frame(spmod))) # works because model.frame(object) works (object has class spmod)
    expect_vector(model.response(model.frame(spmod))) # works because model.frame(object) works (object has class spmod)
  })

  test_that("free generics work autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_true(is.language(terms(spmod))) # works because there is an object spmod$terms
    expect_true(is.call(getCall(spmod))) # works because there is an object spmod$call
    expect_equal(update(formula(spmod), y ~ 1), y ~ 1) # works because formula(object) works (object has class spmod)
    spmod <- update(spmod, y ~ 1, spcov_type = "sar")
    expect_s3_class(spmod, "spautor")
    expect_equal(formula(spmod), y ~ 1)
    expect_s3_class(coefficients(spmod, type = "spcov"), "sar")
    spmod <- update(spmod, . ~ . + offset(x))
    expect_vector(model.offset(model.frame(spmod))) # works because model.frame(object) works (object has class spmod)
    expect_vector(model.response(model.frame(spmod))) # works because model.frame(object) works (object has class spmod)
  })

  ##############################################################################
  ############################ glance (test-glance.R)
  ##############################################################################

  test_that("glance works geostatistical", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_s3_class(glance(spmod), "tbl")
    expect_equal(NROW(glance(spmod)), 1)
    expect_equal(NCOL(glance(spmod)), 9)
    expect_false(any(is.na(glance(spmod))))

    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, estmethod = "sv-wls")
    expect_s3_class(glance(spmod), "tbl")
    expect_equal(NROW(glance(spmod)), 1)
    expect_equal(NCOL(glance(spmod)), 9)
    expect_true(any(is.na(glance(spmod))))
  })

  test_that("glance works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_s3_class(glance(spmod), "tbl")
    expect_equal(NCOL(glance(spmod)), 9)
    expect_false(any(is.na(glance(spmod))))
  })

  ##############################################################################
  ############################ glances (test-glances.R)
  ##############################################################################

  test_that("glances works geostatistical", {
    spmod1 <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    spmod2 <- splm(y ~ x, exdata, "matern", xcoord, ycoord)
    expect_s3_class(glances(spmod1, spmod2), "tbl")
    expect_equal(NROW(glances(spmod1, spmod2)), 2)
    expect_equal(NCOL(glances(spmod1, spmod2)), 10)
    expect_equal(rbind(glance(spmod1), glance(spmod2)), glances(spmod1, spmod2)[, -1], ignore_attr = TRUE)
  })

  test_that("glances works autoregressive", {
    spmod1 <- spautor(y ~ x, exdata_poly, "car")
    spmod2 <- spautor(y ~ x, exdata_poly, "sar")
    expect_s3_class(glances(spmod1, spmod2), "tbl")
    expect_equal(NROW(glances(spmod1, spmod2)), 2)
    expect_equal(NCOL(glances(spmod1, spmod2)), 10)
    expect_equal(rbind(glance(spmod2), glance(spmod1)), glances(spmod1, spmod2)[, -1], ignore_attr = TRUE)
  })

  test_that("glances works iteratively splm spcov_type", {
    spmod <- splm(y ~ x, exdata, c("exponential", "matern"), xcoord, ycoord)
    expect_s3_class(glances(spmod), "tbl")
    expect_equal(NROW(glances(spmod)), 2)
    expect_equal(NCOL(glances(spmod)), 10)
    expect_equal(rbind(glance(spmod$exponential), glance(spmod$matern)), glances(spmod)[, -1], ignore_attr = TRUE)
  })

  test_that("glances works iteratively splm spcov_initial", {
    spcov_init <- lapply(c("exponential", "matern"), function(x) spcov_initial(x, de = 1))
    spmod <- splm(y ~ x, exdata, spcov_initial = spcov_init, xcoord = xcoord, ycoord = ycoord)
    expect_s3_class(glances(spmod), "tbl")
    expect_equal(NROW(glances(spmod)), 2)
    expect_equal(NCOL(glances(spmod)), 10)
    expect_equal(rbind(glance(spmod$spcov_initial_1), glance(spmod$spcov_initial_2)), glances(spmod)[, -1], ignore_attr = TRUE)
  })

  test_that("glances works iteratively spautor spcov_type", {
    spmod <- spautor(y ~ x, exdata_poly, c("car", "sar"))
    expect_s3_class(glances(spmod), "tbl")
    expect_equal(NROW(glances(spmod)), 2)
    expect_equal(NCOL(glances(spmod)), 10)
    expect_equal(rbind(glance(spmod$car), glance(spmod$sar)), glances(spmod, sort_by = "order")[, -1], ignore_attr = TRUE)
    # sort by order here because sort by is off in natural implementation by AICc
  })

  test_that("glances works iteratively spautor spcov_initial", {
    spcov_init <- lapply(c("car", "sar"), function(x) spcov_initial(x, de = 1))
    spmod <- spautor(y ~ x, exdata_poly, spcov_initial = spcov_init)
    expect_s3_class(glances(spmod), "tbl")
    expect_equal(NROW(glances(spmod)), 2)
    expect_equal(NCOL(glances(spmod)), 10)
    expect_equal(rbind(glance(spmod$spcov_initial_1), glance(spmod$spcov_initial_2)), glances(spmod, sort_by = "order")[, -1], ignore_attr = TRUE)
  })

  ##############################################################################
  ############################ hatvalues (test-hatvalues.R)
  ##############################################################################

  test_that("hatvalues works geostatistical", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_vector(hatvalues(spmod))
    expect_true(min(hatvalues(spmod)) >= 0)
    expect_true(max(hatvalues(spmod)) <= 1)
    expect_equal(sum(hatvalues(spmod)), spmod$p)
  })

  test_that("hatvalues works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_vector(hatvalues(spmod))
    expect_true(min(hatvalues(spmod)) >= 0)
    expect_true(max(hatvalues(spmod)) <= 1)
    expect_equal(sum(hatvalues(spmod)), spmod$p)
  })

  ##############################################################################
  ############################ influence (test-influence.R)
  ##############################################################################

  test_that("influence works geostatistical", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_error(influence(spmod), NA)
    expect_equal(NROW(influence(spmod)), 100)
    expect_equal(NCOL(influence(spmod)), 4)
    expect_equal(colnames(influence(spmod)), c(".resid", ".hat", ".cooksd", ".std.resid"))
    expect_s3_class(influence(spmod), "tbl")
  })

  test_that("influence works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_error(influence(spmod), NA)
    expect_equal(NROW(influence(spmod)), 49)
    expect_equal(NCOL(influence(spmod)), 4)
    expect_equal(colnames(influence(spmod)), c(".resid", ".hat", ".cooksd", ".std.resid"))
    expect_s3_class(influence(spmod), "tbl")
  })

  ##############################################################################
  ############################ labels (test-labels.R)
  ##############################################################################

  test_that("labels works geostatistical", {
    spmod <- splm(y ~ x + group, exdata, "exponential", xcoord, ycoord)
    expect_equal(labels(spmod), c("x", "group"))
  })

  test_that("labels works autoregressive", {
    spmod <- spautor(y ~ x + group, exdata_poly, "car")
    expect_equal(labels(spmod), c("x", "group"))
  })

  ##############################################################################
  ############################ logLik (test-logLik.R)
  ##############################################################################

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

  ##############################################################################
  ############################ loocv (test-loocv.R)
  ##############################################################################

  test_that("loocv works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_vector(loocv(spmod))
    expect_vector(loocv(spmod, local = TRUE))
    # cores 2 for cran check
    if (test_local) { ##### local test
      expect_vector(loocv(spmod, local = list(parallel = TRUE, ncores = 2)))
      expect_equal(length(loocv(spmod, cv_predict = TRUE)), 2)
      expect_equal(length(loocv(spmod, cv_predict = TRUE, local = TRUE)), 2)
    }
    expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE)), 3)
    expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = TRUE)), 3)
    if (test_local) { ##### local test
      expect_equal(length(loocv(spmod, se.fit = TRUE)), 2)
      expect_equal(length(loocv(spmod, se.fit = TRUE, local = TRUE)), 2)
    }
    # cores 2 for cran check
    if (test_local) { ##### local test
      expect_equal(length(loocv(spmod, cv_predict = TRUE, local = list(parallel = TRUE, ncores = 2))), 2)
      expect_equal(length(loocv(spmod, cv_predict = TRUE, local = list(parallel = TRUE, method = "all", ncores = 2))), 2)
    }
    expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = list(parallel = TRUE, ncores = 2))), 3)
    expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = list(parallel = TRUE, method = "all", ncores = 2))), 3)
    if (test_local) { ##### local test
      expect_equal(length(loocv(spmod, se.fit = TRUE, local = list(parallel = TRUE, ncores = 2))), 2)
      expect_equal(length(loocv(spmod, se.fit = TRUE, local = list(parallel = TRUE, method = "all", ncores = 2))), 2)
    }

    # random effects
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, random = ~group)
    expect_vector(loocv(spmod))
    expect_vector(loocv(spmod, local = TRUE))
  })

  test_that("loocv works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_vector(loocv(spmod))
    # cores 2 for cran check
    if (test_local) {
      expect_vector(loocv(spmod, local = list(parallel = TRUE, ncores = 2)))
      expect_equal(length(loocv(spmod, cv_predict = TRUE)), 2)
    }
    expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE)), 3)
    # cores 2 for cran check
    if (test_local) {
      expect_equal(length(loocv(spmod, cv_predict = TRUE, local = list(parallel = TRUE, ncores = 2))), 2)
      expect_equal(length(loocv(spmod, se.fit = TRUE, local = list(parallel = TRUE, ncores = 2))), 2)
    }
    expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = list(parallel = TRUE, ncores = 2))), 3)


    # random effects
    spmod <- spautor(y ~ x, exdata_poly, "car", random = ~group)
    expect_vector(loocv(spmod))
    expect_vector(loocv(spmod, local = TRUE))


    # missing data
    spmod <- spautor(y ~ x, exdata_Mpoly, "car")
    expect_vector(loocv(spmod))
    # cores 2 for cran check
    if (test_local) {
      expect_vector(loocv(spmod, local = list(parallel = TRUE, ncores = 2)))
      expect_equal(length(loocv(spmod, cv_predict = TRUE)), 2)
    }
    expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE)), 3)
    # cores 2 for cran check
    if (test_local) {
      expect_equal(length(loocv(spmod, cv_predict = TRUE, local = list(parallel = TRUE, ncores = 2))), 2)
      expect_equal(length(loocv(spmod, se.fit = TRUE, local = list(parallel = TRUE, ncores = 2))), 2)
    }
    expect_equal(length(loocv(spmod, cv_predict = TRUE, se.fit = TRUE, local = list(parallel = TRUE, ncores = 2))), 3)

  })

  ##############################################################################
  ############################ model.frame (test-model.frame.R)
  ##############################################################################

  test_that("model.frame works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_equal(NCOL(model.frame(spmod)), 2)
    expect_equal(NROW(model.frame(spmod)), 100)
    expect_equal(model.frame(spmod), model.frame(y ~ x, exdata, drop.unused.levels = TRUE, na.action = na.omit))
  })

  test_that("model.frame works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_equal(NCOL(model.frame(spmod)), 2)
    expect_equal(NROW(model.frame(spmod)), 49)
    expect_equal(model.frame(spmod), model.frame(y ~ x, exdata_poly, drop.unused.levels = TRUE, na.action = na.omit))
  })

  ##############################################################################
  ############################ model.matrix (test-model.matrix.R)
  ##############################################################################

  test_that("model.matrix works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_equal(NCOL(model.matrix(spmod)), 2)
    expect_equal(NROW(model.matrix(spmod)), 100)
    expect_equal(model.matrix(spmod), model.matrix(y ~ x, model.frame(y ~ x, exdata, drop.unused.levels = TRUE, na.action = na.omit)))
  })

  test_that("model.matrix works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_equal(NCOL(model.matrix(spmod)), 2)
    expect_equal(NROW(model.matrix(spmod)), 49)
    expect_equal(model.matrix(spmod), model.matrix(y ~ x, model.frame(y ~ x, exdata_poly, drop.unused.levels = TRUE, na.action = na.omit)))
  })

  ##############################################################################
  ############################ plot (test-plot.R)
  ##############################################################################

  test_that("plot works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)

    # plot 1
    expect_error(plot(spmod, which = 1), NA)

    # plot 2
    expect_error(plot(spmod, which = 2), NA)

    # plot 3
    expect_error(plot(spmod, which = 3), NA)

    # plot 4
    expect_error(plot(spmod, which = 4), NA)

    # plot 5
    expect_error(plot(spmod, which = 5), NA)

    # plot 6
    expect_error(plot(spmod, which = 6), NA)

    # plot 7
    expect_error(plot(spmod, which = 7), NA)

    # plot 8
    expect_error(plot(spmod, which = 8), NA)

    # plot 9 (return error)
    expect_error(plot(spmod, which = 9))

    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, anisotropy = TRUE)

    # plot 8
    expect_error(plot(spmod, which = 8), NA)
  })

  test_that("plot works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")

    # plot 1
    expect_error(plot(spmod, which = 1), NA)

    # plot 2
    expect_error(plot(spmod, which = 2), NA)

    # plot 3
    expect_error(plot(spmod, which = 3), NA)

    # plot 4
    expect_error(plot(spmod, which = 4), NA)

    # plot 5
    expect_error(plot(spmod, which = 5), NA)

    # plot 6
    expect_error(plot(spmod, which = 6), NA)

    # plot 7
    expect_error(plot(spmod, which = 7))

    # plot 8
    expect_error(plot(spmod, which = 8))
  })

  ##############################################################################
  ############################ predict (test-predict.R)
  ##############################################################################

  test_that("Prediction for splm works", {
    spcov_type <- "exponential"
    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
    expect_error(predict(smod, newexdata), NA)
    expect_error(predict(smod, newexdata, interval = "prediction"), NA)
    expect_error(predict(smod, newexdata, interval = "confidence"), NA)
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))
    expect_true(all(predict(smod, newexdata, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction for splm works anisotropy", {
    spcov_type <- "exponential"
    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", anisotropy = TRUE)
    expect_error(predict(smod, newexdata), NA)
    expect_error(predict(smod, newexdata, interval = "prediction"), NA)
    expect_error(predict(smod, newexdata, interval = "confidence"), NA)
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))
    expect_true(all(predict(smod, newexdata, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction for splm works with random effects", {
    spcov_type <- "exponential"
    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group)
    expect_error(predict(smod, newexdata), NA)
    expect_error(predict(smod, newexdata, interval = "prediction"), NA)
    expect_error(predict(smod, newexdata, interval = "confidence"), NA)
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))
    expect_true(all(predict(smod, newexdata, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction for splm works with partition factor", {
    spcov_type <- "exponential"
    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", partition_factor = ~group)
    expect_error(predict(smod, newexdata), NA)
    expect_error(predict(smod, newexdata, interval = "prediction"), NA)
    expect_error(predict(smod, newexdata, interval = "confidence"), NA)
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))
    expect_true(all(predict(smod, newexdata, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction works for big data", {
    spcov_type <- "exponential"
    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
    expect_error(predict(smod, newexdata), NA)
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))
    # big data
    expect_error(predict(smod, newexdata, local = TRUE), NA)
    expect_error(predict(smod, newexdata, interval = "prediction", local = TRUE), NA)
    expect_error(predict(smod, newexdata, interval = "confidence", local = TRUE), NA)
    # CRAN FIXES CORES AT 2 MAX
    expect_error(predict(smod, newexdata, local = list(parallel = TRUE, ncores = 2)), NA)
    expect_error(predict(smod, newexdata, interval = "prediction", local = list(parallel = TRUE, ncores = 2)), NA)
    expect_error(predict(smod, newexdata, interval = "confidence", local = list(parallel = TRUE, ncores = 2)), NA)
    expect_equal(length(predict(smod, newexdata, local = TRUE)), NROW(newexdata))
    expect_true(all(predict(smod, newexdata, local = TRUE, se.fit = TRUE)$se.fit >= 0))
    expect_error(predict(smod, newexdata, local = list(method = "distance")), NA)
    expect_equal(length(predict(smod, newexdata, local = list(method = "distance"))), NROW(newexdata))
    expect_error(predict(smod, newexdata, local = list(method = "covariance")), NA)
    expect_equal(length(predict(smod, newexdata, local = list(method = "covariance"))), NROW(newexdata))
    expect_error(predict(smod, newexdata, local = list(method = "distance", size = 10)), NA)
    expect_equal(length(predict(smod, newexdata, local = list(method = "distance", size = 10))), NROW(newexdata))
    expect_error(predict(smod, newexdata, local = list(method = "covariance", size = 10)), NA)
    expect_equal(length(predict(smod, newexdata, local = list(method = "covariance", size = 10))), NROW(newexdata))
  })

  test_that("Prediction for splm works for missing data", {
    spcov_type <- "exponential"
    smod <- splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
    expect_error(predict(smod), NA)
    expect_error(predict(smod, interval = "prediction"), NA)
    expect_error(predict(smod, interval = "confidence"), NA)
    expect_equal(length(predict(smod)), sum(is.na(exdata_M$y)))
    expect_true(all(predict(smod, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction for spautor works", {
    spcov_type <- "car"
    smod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "reml")
    expect_error(predict(smod), NA)
    expect_error(predict(smod, interval = "prediction"), NA)
    expect_error(predict(smod, interval = "confidence"), NA)
    expect_equal(length(predict(smod)), sum(is.na(exdata_Mpoly$y)))
    expect_true(all(predict(smod, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction for spautor works parallel", {
    spcov_type <- "car"
    smod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "reml")
    # CRAN FIXES CORES AT 2 MAX
    expect_error(predict(smod, local = list(parallel = TRUE, ncores = 2)), NA)
    expect_error(predict(smod, interval = "prediction", local = list(parallel = TRUE, ncores = 2)), NA)
    expect_error(predict(smod, interval = "confidence", local = list(parallel = TRUE, ncores = 2)), NA)
    expect_equal(length(predict(smod)), sum(is.na(exdata_Mpoly$y)))
    expect_true(all(predict(smod, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction for spautor works with random effects", {
    spcov_type <- "car"
    smod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "reml", random = ~group)
    expect_error(predict(smod), NA)
    expect_error(predict(smod, interval = "prediction"), NA)
    expect_error(predict(smod, interval = "confidence"), NA)
    expect_equal(length(predict(smod)), sum(is.na(exdata_Mpoly$y)))
    expect_true(all(predict(smod, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction for spautor works with partition factor", {
    spcov_type <- "car"
    smod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "reml", partition = ~group)
    expect_error(predict(smod), NA)
    expect_error(predict(smod, interval = "prediction"), NA)
    expect_error(predict(smod, interval = "confidence"), NA)
    expect_equal(length(predict(smod)), sum(is.na(exdata_Mpoly$y)))
    expect_true(all(predict(smod, se.fit = TRUE)$se.fit >= 0))
  })

  test_that("Prediction works for other covariances", {
    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "gaussian", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, spcov_type = "triangular", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "circular", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "cubic", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "pentaspherical", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, spcov_type = "cosine", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "wave", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "jbessel", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "gravity", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "rquad", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "magnetic", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "cauchy", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "pexponential", estmethod = "reml")
    expect_equal(length(predict(smod, newexdata)), NROW(newexdata))

    smod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car", estmethod = "reml")
    expect_equal(length(predict(smod)), sum(is.na(exdata_Mpoly$y)))

    smod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "sar", estmethod = "reml")
    expect_equal(length(predict(smod)), sum(is.na(exdata_Mpoly$y)))
  })

  test_that("errors occur", {
    spcov_type <- "exponential"
    spmod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
    expect_error(predict(spmod))
    expect_error(predict(spmod, newexdata = newexdata, local = list(method = "xyz")))

    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_error(predict(spmod))
  })

  test_that("prediction values match for both approaches", {
    spmod1 <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    pred1 <- predict(spmod1, newdata = newexdata)
    newexdata$y <- NA
    exdata_with_NA <- rbind(exdata, newexdata)
    spmod2 <- splm(y ~ x, exdata_with_NA, "exponential", xcoord, ycoord)
    pred2 <- predict(spmod2)
    pred3 <- predict(spmod2, newdata = spmod2$newdata)

    spmod1$call <- NULL # calls are different among two splm() calls
    spmod2$call <- NULL
    names(pred1) <- NULL # names start at 1
    names(pred2) <- NULL # names start at index in data
    names(pred3) <- NULL # names start at 1
    expect_equal(summary(spmod1), summary(spmod2))
    expect_equal(pred1, pred2)
    expect_equal(pred2, pred3)


    spmod1 <- splm(y ~ poly(x, degree = 2, raw = TRUE), exdata, "exponential", xcoord, ycoord)
    pred1 <- predict(spmod1, newdata = newexdata)
    newexdata$y <- NA
    exdata_with_NA <- rbind(exdata, newexdata)
    spmod2 <- splm(y ~ poly(x, degree = 2, raw = TRUE), exdata_with_NA, "exponential", xcoord, ycoord)
    pred2 <- predict(spmod2)
    pred3 <- predict(spmod2, newdata = spmod2$newdata)

    spmod1$call <- NULL # calls are different among two splm() calls
    spmod2$call <- NULL
    names(pred1) <- NULL # names start at 1
    names(pred2) <- NULL # names start at index in data
    names(pred3) <- NULL # names start at 1
    expect_equal(summary(spmod1), summary(spmod2))
    expect_equal(pred1, pred2)
    expect_equal(pred2, pred3)

    spmod1 <- splm(y ~ poly(x, degree = 2, raw = FALSE), exdata, "exponential", xcoord, ycoord)
    pred1 <- predict(spmod1, newdata = newexdata)
    newexdata$y <- NA
    exdata_with_NA <- rbind(exdata, newexdata)
    spmod2 <- splm(y ~ poly(x, degree = 2, raw = FALSE), exdata_with_NA, "exponential", xcoord, ycoord)
    pred2 <- predict(spmod2)
    pred3 <- predict(spmod2, newdata = spmod2$newdata)

    spmod1$call <- NULL # calls are different among two splm() calls
    spmod2$call <- NULL
    names(pred1) <- NULL # names start at 1
    names(pred2) <- NULL # names start at index in data
    names(pred3) <- NULL # names start at 1
    expect_equal(summary(spmod1), summary(spmod2))
    expect_equal(pred1, pred2)
    expect_equal(pred2, pred3)
  })

  test_that("prediction values match for both and lm comparison", {

    # no poly
    spmod1 <- splm(y ~ x, exdata, "none")
    pred1 <- predict(spmod1, newdata = newexdata)
    newexdata$y <- NA
    exdata_with_NA <- rbind(exdata, newexdata)
    spmod2 <- splm(y ~ x, exdata_with_NA, "none")
    pred2 <- predict(spmod2)
    pred3 <- predict(spmod2, newdata = spmod2$newdata)

    spmod1$call <- NULL # calls are different among two splm() calls
    spmod2$call <- NULL
    names(pred1) <- NULL # names start at 1
    names(pred2) <- NULL # names start at index in data
    names(pred3) <- NULL # names start at 1
    expect_equal(summary(spmod1), summary(spmod2))
    expect_equal(pred1, pred2)
    expect_equal(pred2, pred3)

    ## compare lm
    lmod1 <- lm(y ~ x, exdata)
    lmpred1 <- predict(lmod1, newexdata)
    expect_equal(unname(pred1), unname(lmpred1))

    # poly raw
    spmod1 <- splm(y ~ poly(x, degree = 2, raw = TRUE), exdata, "none")
    pred1 <- predict(spmod1, newdata = newexdata)
    newexdata$y <- NA
    exdata_with_NA <- rbind(exdata, newexdata)
    spmod2 <- splm(y ~ poly(x, degree = 2, raw = TRUE), exdata_with_NA, "none")
    pred2 <- predict(spmod2)
    pred3 <- predict(spmod2, newdata = spmod2$newdata)

    spmod1$call <- NULL # calls are different among two splm() calls
    spmod2$call <- NULL
    names(pred1) <- NULL # names start at 1
    names(pred2) <- NULL # names start at index in data
    names(pred3) <- NULL # names start at 1
    expect_equal(summary(spmod1), summary(spmod2))
    expect_equal(pred1, pred2)
    expect_equal(pred2, pred3)

    ## compare lm
    lmod1 <- lm(y ~ poly(x, degree = 2, raw = TRUE), exdata)
    lmpred1 <- predict(lmod1, newexdata)
    expect_equal(unname(pred1), unname(lmpred1))

    # poly no raw
    spmod1 <- splm(y ~ poly(x, degree = 2, raw = FALSE), exdata, "none")
    pred1 <- predict(spmod1, newdata = newexdata)
    newexdata$y <- NA
    exdata_with_NA <- rbind(exdata, newexdata)
    spmod2 <- splm(y ~ poly(x, degree = 2, raw = FALSE), exdata_with_NA, "none")
    pred2 <- predict(spmod2)
    pred3 <- predict(spmod2, newdata = spmod2$newdata)

    spmod1$call <- NULL # calls are different among two splm() calls
    spmod2$call <- NULL
    names(pred1) <- NULL # names start at 1
    names(pred2) <- NULL # names start at index in data
    names(pred3) <- NULL # names start at 1
    expect_equal(summary(spmod1), summary(spmod2))
    expect_equal(pred1, pred2)
    expect_equal(pred2, pred3)

    ## compare lm
    lmod1 <- lm(y ~ poly(x, degree = 2, raw = FALSE), exdata)
    lmpred1 <- predict(lmod1, newexdata)
    expect_equal(unname(pred1), unname(lmpred1))
  })

  test_that("prediction values match for both approaches autoregressive", {
    spmod1 <- spautor(y ~ poly(x, degree = 2, raw = TRUE), exdata_Mpoly, "car")
    expect_error(predict(spmod1), NA)

    spmod1 <- spautor(y ~ poly(x, degree = 2, raw = FALSE), exdata_Mpoly, "car")
    expect_error(predict(spmod1), NA)
  })

  test_that("prediction no error with order 2 polynomial one prediction row", {
    # there is a bug in lm() trying to do the same thing
    spmod <- splm(y ~ poly(xcoord, ycoord, degree = 1), exdata, "none")
    expect_error(predict(spmod, newdata = newexdata), NA)
    expect_error(predict(spmod, newdata = newexdata[1, , drop = FALSE]), NA)
    expect_equal(predict(spmod, newdata = newexdata)[[1]], predict(spmod, newdata = newexdata[1, , drop = FALSE])[[1]])

    # there is a bug in lm() trying to do the same thing
    spmod <- splm(y ~ poly(xcoord, ycoord, degree = 1), exdata, "exponential", xcoord, ycoord)
    expect_error(predict(spmod, newdata = newexdata), NA)
    expect_error(predict(spmod, newdata = newexdata[1, , drop = FALSE]), NA)
    expect_equal(predict(spmod, newdata = newexdata)[[1]], predict(spmod, newdata = newexdata[1, , drop = FALSE])[[1]])

    # there is a bug in lm() trying to do the same thing
    exdata_Mpoly$x2 <- rnorm(NROW(exdata_Mpoly))
    spmod <- spautor(y ~ poly(x, x2, degree = 1), exdata_Mpoly, "car")
    expect_error(predict(spmod), NA)
  })

  test_that("prediction works when newdata does not have all factor levels", {
    # geostatistical
    spmod <- splm(y ~ group, exdata, "exponential", xcoord, ycoord)
    newexdata_sub <- newexdata[1, , drop = FALSE]
    newexdata_sub$group <- as.character(newexdata_sub$group)
    expect_error(predict(spmod, newexdata_sub), NA)

    # autoregressive
    exdata_Mpoly_new <- exdata_Mpoly
    exdata_Mpoly_new$group <- as.character(exdata_Mpoly_new$group)
    spmod <- spautor(y ~ group, exdata_Mpoly_new, "car")
    expect_error(predict(spmod), NA)
  })

  test_that("prediction for splmRF works", {
    spcov_type <- "exponential"

    sprfmod <- splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
    expect_vector(predict(sprfmod, newdata = newexdata))

    exdata$y[1] <- NA
    sprfmod <- splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
    expect_vector(predict(sprfmod))


    spcov_type <- c("exponential", "matern")

    sprfmod <- splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
    expect_vector(predict(sprfmod, newdata = newexdata)[[1]])
    expect_vector(predict(sprfmod, newdata = newexdata)[[2]])

    exdata$y[1] <- NA
    sprfmod <- splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
    expect_vector(predict(sprfmod)[[1]])
    expect_vector(predict(sprfmod)[[2]])
  })

  test_that("prediction for spautorRF works", {
    spcov_type <- "car"
    sprfmod <- spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type)
    expect_vector(predict(sprfmod))

    spcov_type <- c("car", "sar")
    sprfmod <- spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type)
    expect_vector(predict(sprfmod))
    expect_vector(predict(sprfmod)[[1]])
    expect_vector(predict(sprfmod)[[2]])
  })

  test_that("prediction for splm offset works", {
    exdata$offset <- 2
    exdata$y2 <- exdata$y - exdata$offset
    newexdata$offset <- 2
    spmod1 <- splm(y ~ x + offset(offset), exdata, "exponential", xcoord, ycoord)
    spmod2 <- splm(y2 ~ x,  exdata, "exponential", xcoord, ycoord)
    expect_equal(fitted(spmod1), fitted(spmod2) + exdata$offset)
  })

  test_that("prediction for spautor offset works", {
    exdata_Mpoly$offset <- 2
    exdata_Mpoly$y2 <- exdata_Mpoly$y - exdata_Mpoly$offset
    spmod1 <- spautor(y ~ x + offset(offset), exdata_Mpoly, "car")
    spmod2 <- spautor(y2 ~ x,  exdata_Mpoly, "car")
    expect_equal(fitted(spmod1), fitted(spmod2) + exdata_Mpoly$offset[!is.na(exdata_Mpoly$y)])
  })

  ##############################################################################
  ############################ print (test-print.R)
  ##############################################################################

  test_that("printing works for geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_output(print(spmod))
    expect_output(print(summary(spmod)))
    expect_output(print(anova(spmod)))
  })

  test_that("printing works for auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_output(print(spmod))
    expect_output(print(summary(spmod)))
    expect_output(print(anova(spmod)))
  })

  ##############################################################################
  ############################ pseudoR2 (test-pseudoR2.R)
  ##############################################################################

  test_that("pseudoR2 works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_vector(pseudoR2(spmod))
    expect_vector(pseudoR2(spmod, adjust = TRUE))
  })

  test_that("pseudoR2 works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_vector(pseudoR2(spmod))
    expect_vector(pseudoR2(spmod, adjust = TRUE))
  })

  test_that("pseudoR2 matches classical r-squared", {
    lmod <- lm(y ~ x, exdata)
    spmod <- splm(y ~ x, exdata, "none")
    expect_equal(summary(lmod)$r.squared, pseudoR2(spmod))
    expect_equal(summary(lmod)$adj.r.squared, pseudoR2(spmod, adjust = TRUE))
  })

  ##############################################################################
  ############################ residuals (test-residuals.R)
  ##############################################################################

  test_that("residuals works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_vector(residuals(spmod))
    expect_equal(residuals(spmod), exdata$y - fitted(spmod))
    expect_vector(residuals(spmod, type = "response"))
    expect_vector(residuals(spmod, type = "pearson"))
    expect_vector(residuals(spmod, type = "standardized"))
    expect_equal(residuals(spmod, type = "standardized"), rstandard(spmod))
  })

  test_that("residuals works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_vector(residuals(spmod))
    expect_equal(residuals(spmod), exdata_poly$y - fitted(spmod))
    expect_vector(residuals(spmod, type = "response"))
    expect_vector(residuals(spmod, type = "pearson"))
    expect_vector(residuals(spmod, type = "standardized"))
    expect_equal(residuals(spmod, type = "standardized"), rstandard(spmod))
  })

  ################################################################################
  ############################ sp_objects (test-sp_objects.R)
  ################################################################################

  test_that("sp object error tests geostatistical", {
    # first we make "fake" sp objects that can be loaded without also loading sp
    exdata_sp <- exdata
    attr(class(exdata_sp), "package") <- "sp"
    expect_error(splm(y ~ x, exdata_sp, "exponential", xcoord = xcoord, ycoord = ycoord))
    expect_error(esv(y ~ x, exdata_sp, xcoord = xcoord, ycoord = ycoord))
    spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
    expect_error(sprnorm(spcov_params_val, data = exdata_sp, xcoord = xcoord, ycoord = ycoord))

    spmod <- splm(y ~ x, exdata, "exponential", xcoord = xcoord, ycoord = ycoord)
    newexdata_sp <- newexdata
    attr(class(newexdata_sp), "package") <- "sp"
    expect_error(predict(spmod, newexdata_sp))
    expect_error(augment(spmod, newdata = newexdata_sp))
  })

  test_that("sp object error tests autoregressive", {
    # first we make "fake" sp objects that can be loaded without also loading sp
    exdata_poly_sp <- exdata_poly
    attr(class(exdata_poly_sp), "package") <- "sp"
    expect_error(spautor(y ~ x, exdata_poly_sp, "car"))
    spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 0.5, extra = 0)
    expect_error(sprnorm(spcov_params_val, data = exdata_poly_sp))
  })

  ##############################################################################
  ############################ sprnorm (test-sprnorm.R)
  ##############################################################################

  test_that("the simulation runs for exponential", {
    spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for exponential anisotropy", {
    spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord, anisotropy = TRUE))
  })

  test_that("the simulation runs for exponential random effect", {
    spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
    randcov_params_val <- randcov_params(group = 1)
    expect_vector(sprnorm(spcov_params_val,
                          data = exdata, xcoord = xcoord, ycoord = ycoord,
                          random = ~group, randcov_params = randcov_params_val
    ))
  })

  test_that("the simulation runs for exponential partitioning", {
    spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord, partition_factor = ~group))
  })

  test_that("the simulation runs for exponential (matrix)", {
    spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
    sprnorm_matrix <- sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord, samples = 5)
    expect_equal(nrow(sprnorm_matrix), nrow(exdata))
    expect_equal(ncol(sprnorm_matrix), 5)
  })

  test_that("the simulation runs for car", {
    spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 1 / 2, extra = 0)
    expect_vector(sprnorm(spcov_params_val, data = exdata_poly))
  })

  test_that("the simulation runs for car unconnected", {
    spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 1 / 2, extra = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata_Upoly))
  })

  test_that("the simulation runs for none random effect", {
    spcov_params_val <- spcov_params("none", ie = 1)
    randcov_params_val <- randcov_params(group = 1)
    expect_vector(sprnorm(spcov_params_val,
                          data = exdata, xcoord = xcoord, ycoord = ycoord,
                          random = ~group, randcov_params = randcov_params_val
    ))
  })

  test_that("the simulation runs for car random effects", {
    spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 1 / 2, extra = 0)
    expect_vector(sprnorm(spcov_params_val, data = exdata_poly))
    randcov_params_val <- randcov_params(group = 1)
    expect_vector(sprnorm(spcov_params_val,
                          data = exdata_poly, xcoord = xcoord, ycoord = ycoord,
                          random = ~group, randcov_params = randcov_params_val
    ))
  })

  test_that("the simulation runs for car partition factor", {
    spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 1 / 2, extra = 0)
    expect_vector(sprnorm(spcov_params_val, data = exdata_poly, partition_factor = ~group))
  })




  test_that("the simulation runs for spherical", {
    spcov_params_val <- spcov_params("spherical", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for gaussian", {
    spcov_params_val <- spcov_params("gaussian", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for triangular", {
    spcov_params_val <- spcov_params("triangular", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord))
  })

  test_that("the simulation runs for circular", {
    spcov_params_val <- spcov_params("circular", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for cubic", {
    spcov_params_val <- spcov_params("cubic", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for pentaspherical", {
    spcov_params_val <- spcov_params("pentaspherical", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for cosine", {
    spcov_params_val <- spcov_params("cosine", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord))
  })

  test_that("the simulation runs for wave", {
    spcov_params_val <- spcov_params("wave", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for jbessel", {
    spcov_params_val <- spcov_params("jbessel", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for gravity", {
    spcov_params_val <- spcov_params("gravity", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for rquad", {
    spcov_params_val <- spcov_params("rquad", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for magnetic", {
    spcov_params_val <- spcov_params("magnetic", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for none", {
    spcov_params_val <- spcov_params("none", ie = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4))
  })

  test_that("the simulation runs for matern", {
    spcov_params_val <- spcov_params("matern", de = 1, ie = 1, range = 1, extra = 1 / 2)
    expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for cauchy", {
    spcov_params_val <- spcov_params("cauchy", de = 1, ie = 1, range = 1, extra = 1 / 2)
    expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for pexponential", {
    spcov_params_val <- spcov_params("pexponential", de = 1, ie = 1, range = 1, extra = 1 / 2)
    expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  })

  test_that("the simulation runs for sar", {
    spcov_params_val <- spcov_params("sar", de = 1, ie = 1, range = 1 / 2, extra = 0)
    expect_vector(sprnorm(spcov_params_val, data = exdata_poly))
  })

  test_that("the simulation runs for sar unconnected", {
    spcov_params_val <- spcov_params("sar", de = 1, ie = 1, range = 1 / 2, extra = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata_Upoly))
  })

  test_that("quoting works", {
    spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = "xcoord", ycoord = "ycoord"))
  })

  test_that("sf works", {
    exdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"))
    spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
    expect_vector(sprnorm(spcov_params_val, data = exdata_sf))
  })


  ##############################################################################
  ############################ summary (test-summary.R)
  ##############################################################################

  test_that("summary works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_error(summary(spmod), NA)
  })

  test_that("summary works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_error(summary(spmod), NA)
  })

  ##############################################################################
  ############################ tidy (test-tidy.R)
  ##############################################################################

  test_that("tidy works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_s3_class(tidy(spmod), "tbl")
    expect_equal(ncol(tidy(spmod)), 5)
    expect_s3_class(tidy(spmod, effects = "spcov"), "tbl")
    expect_equal(ncol(tidy(spmod, effects = "spcov")), 3)
    expect_null(tidy(spmod, effects = "randcov"))
    expect_s3_class(tidy(anova(spmod)), "tbl")
    expect_equal(ncol(tidy(anova(spmod))), 4)
  })

  test_that("tidy works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, random = ~group)
    expect_s3_class(tidy(spmod), "tbl")
    expect_equal(ncol(tidy(spmod)), 5)
    expect_s3_class(tidy(spmod, effects = "spcov"), "tbl")
    expect_equal(ncol(tidy(spmod, effects = "spcov")), 3)
    expect_s3_class(tidy(spmod, effects = "randcov"), "tbl")
    expect_equal(ncol(tidy(spmod, effects = "randcov")), 3)
    expect_s3_class(tidy(anova(spmod)), "tbl")
    expect_equal(ncol(tidy(anova(spmod))), 4)
  })

  test_that("tidy works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_s3_class(tidy(spmod), "tbl")
    expect_equal(ncol(tidy(spmod)), 5)
    expect_s3_class(tidy(spmod, effects = "spcov"), "tbl")
    expect_equal(ncol(tidy(spmod, effects = "spcov")), 3)
    expect_null(tidy(spmod, effects = "randcov"))
    expect_s3_class(tidy(anova(spmod)), "tbl")
    expect_equal(ncol(tidy(anova(spmod))), 4)
  })

  test_that("tidy works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car", random = ~group)
    expect_s3_class(tidy(spmod), "tbl")
    expect_equal(ncol(tidy(spmod)), 5)
    expect_s3_class(tidy(spmod, effects = "spcov"), "tbl")
    expect_equal(ncol(tidy(spmod, effects = "spcov")), 3)
    expect_s3_class(tidy(spmod, effects = "randcov"), "tbl")
    expect_equal(ncol(tidy(spmod, effects = "randcov")), 3)
    expect_s3_class(tidy(anova(spmod)), "tbl")
    expect_equal(ncol(tidy(anova(spmod))), 4)
  })

  ##############################################################################
  ############################ varcomp (test-varcomp.R)
  ##############################################################################

  test_that("varcomp works geostatistical", {
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val), 3)
    expect_equal(sum(varcomp_val$proportion), 1)
  })

  test_that("varcomp works geostatistical random", {
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord,
                  ycoord = ycoord, estmethod = "reml", random = ~ group)
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val), 4)
    expect_equal(sum(varcomp_val$proportion), 1)
  })

  test_that("varcomp works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val), 3)
    expect_equal(sum(varcomp_val$proportion), 1)
  })

  test_that("varcomp works autoregressive random", {
    spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml", random = ~ group)
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val), 4)
    expect_equal(sum(varcomp_val$proportion), 1)
  })

  test_that("varcomp works unconnected autoregressive", {
    spmod <- spautor(y ~ x, exdata_Upoly, spcov_type = "car", estmethod = "reml")
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val$connected), 3)
    expect_equal(sum(varcomp_val$connected$proportion), 1)
    expect_equal(NROW(varcomp_val$unconnected), 2)
    expect_equal(sum(varcomp_val$unconnected$proportion), 1)
    expect_equal(length(varcomp_val$ratio), 1)
  })

  test_that("varcomp works unconnected autoregressive random", {
    spmod <- spautor(y ~ x, exdata_Upoly, spcov_type = "car", estmethod = "reml", random = ~ group)
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val$connected), 4)
    expect_equal(sum(varcomp_val$connected$proportion), 1)
    expect_equal(NROW(varcomp_val$unconnected), 3)
    expect_equal(sum(varcomp_val$unconnected$proportion), 1)
    expect_equal(length(varcomp_val$ratio), 1)
  })

  ##############################################################################
  ############################ vcov (test-vcov.R)
  ##############################################################################

  test_that("vcov works geo", {
    spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
    expect_true(isSymmetric(vcov(spmod)))
    expect_true(all(diag(vcov(spmod)) >= 0))
  })

  test_that("vcov works auto", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_true(isSymmetric(vcov(spmod)))
    expect_true(all(diag(vcov(spmod)) >= 0))
  })
}
