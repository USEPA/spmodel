test_that("blank test", {
  expect_null(NULL)
})

test_local <- FALSE # FALSE for CRAN

if (test_local) {
  test_that("exact output match", {

    load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
    load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3943, -0.0730))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.6361, -0.9631))

    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
                  anisotropy = TRUE, random = ~ group, partition_factor = ~ group)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2835, -0.0622))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.2311, -0.2132))

    spmod <- spglm(abs(y) ~ x, exdata, family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3349, 0.0437))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.3043, -0.3981))

    spmod <- spglm(abs(y) ~ x, exdata, family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
                  anisotropy = TRUE, random = ~ group, partition_factor = ~ group)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3349, 0.0437))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.3043, -0.3981)) # estimates same as above models to digit 4 but se differ


    load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

    spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car")
    preds <- predict(spmod, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1543, -0.1212))
    expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2868))

    spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car", random = ~ group, partition_factor = ~ group)
    preds <- predict(spmod, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1481, -0.1226))
    expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2807))

    spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly, spcov_type = "car")
    preds <- predict(spmod, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2769, 0.0317))
    expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2277))

    spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly, spcov_type = "car", random = ~ group, partition_factor = ~ group)
    preds <- predict(spmod, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2781, 0.0305))
    expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2106))

    set.seed(0)
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = list(size = 20))
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE, local = list(size = 20))
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2571, -0.0759))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.5774, -0.9206))

    set.seed(1)
    spmod <- spglm(abs(y) ~ x, exdata, family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = list(size = 20))
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE, local = list(size = 20))
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3556, 0.0501))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.2969, -0.4028))
  })
}

