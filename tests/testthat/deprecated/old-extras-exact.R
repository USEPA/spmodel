test_that("blank test", {
  expect_null(NULL)
})

test_local <- FALSE # FALSE for CRAN

if (test_local) {
  test_that("exact output match", {

    load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
    load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
    load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

    # linear models

    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3943, -0.0730))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.6361, -0.9631))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.3839, 0.0736))
    expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.7647, 0.8270))

    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
                  anisotropy = TRUE, random = ~ group, partition_factor = ~ group)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2835, -0.0622))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.2311, -0.2132))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1111, 0.0826))
    expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.8721, 0.9086))

    spmod <- spglm(abs(y) ~ x, exdata, family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3349, 0.0437))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.3043, -0.3981))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.2744, 0.0853))
    expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.1059, 0.1591))

    spmod <- spglm(abs(y) ~ x, exdata, family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
                  anisotropy = TRUE, random = ~ group, partition_factor = ~ group)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3349, 0.0437))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.3043, -0.3981))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.0918, 0.0853))
    expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.1064, 0.1594))

    # autoregressive models

    spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car")
    preds <- predict(spmod, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1543, -0.1212))
    expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2868))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1532, 0.1269))
    expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(1.0137))

    spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car", random = ~ group, partition_factor = ~ group)
    preds <- predict(spmod, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1481, -0.1226))
    expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2807))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1442, 0.1275))
    expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(1.0147))

    spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly, spcov_type = "car")
    preds <- predict(spmod, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2769, 0.0317))
    expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2277))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1119, 0.1049))
    expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(0.1932))

    spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly, spcov_type = "car", random = ~ group, partition_factor = ~ group)
    preds <- predict(spmod, se.fit = TRUE)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2781, 0.0305))
    expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2106))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1224, 0.1046))
    expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(0.1968))

    ##############
    #### local
    ##############
    local_list <- list(size = 20)
    seed <- 0
    n_skip <- 1

    # linear models

    set.seed(seed)
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = local_list)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE, local = local_list)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2571, -0.0759))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.5774, -0.9206))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.3640, 0.0739))
    expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.7603, 0.8037))
    seed <- seed + n_skip

    set.seed(seed)
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
                  anisotropy = TRUE, random = ~ group, partition_factor = ~ group, local = local_list)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE, local = local_list)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2735, 0.0030))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.0944, -0.3359))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1236, 0.0687))
    expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.8732, 0.9020))
    seed <- seed + n_skip

    set.seed(seed)
    spmod <- spglm(abs(y) ~ x, exdata, family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = local_list)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE, local = local_list)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3556, 0.0501))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.2969, -0.4028))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1255, 0.0851))
    expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.2799, 0.3178))
    seed <- seed + n_skip

    set.seed(seed)
    spmod <- spglm(abs(y) ~ x, exdata, family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
                   anisotropy = TRUE, random = ~ group, partition_factor = ~ group, local = local_list)
    preds <- predict(spmod, newdata = newexdata, se.fit = TRUE, local = local_list)
    expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3349, 0.0437))
    expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.3043, -0.3981))
    expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.0938, 0.0852))
    expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.2874, 0.2940))

  })
}

