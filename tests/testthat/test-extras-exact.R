skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

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

  # optim() does not converge for this fixture (anisotropy + random +
  # partition_factor is numerically marginal); coefficients are still
  # deterministic
  spmod <- suppressWarnings(splm(y ~ x, exdata,
    spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    anisotropy = TRUE, random = ~group, partition_factor = ~group
  ))
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2835, -0.0622))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.2311, -0.2132))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1111, 0.0826))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.8721, 0.9086))

  # a different covariance family than the exponential cases above, to widen
  # the branch coverage this regression suite exercises
  spmod <- splm(y ~ x, exdata, spcov_type = "matern", xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3728, -0.0796))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.6104, -1.1435))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.3440, 0.0723))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.7535, 0.8017))

  spmod <- spglm(abs(y) ~ x, exdata, family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3349, 0.0437))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.3043, -0.3981))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.2744, 0.0853))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.1064, 0.1594))

  spmod <- spglm(abs(y) ~ x, exdata,
    family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    anisotropy = TRUE, random = ~group, partition_factor = ~group
  )
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3349, 0.0437))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.3043, -0.3981))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.0918, 0.0853))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.1064, 0.1595))

  # a different GLM family than the Gamma cases above (a count response)
  exdata_pois <- exdata
  exdata_pois$count <- round(abs(exdata_pois$y) * 3)
  newexdata_pois <- newexdata
  newexdata_pois$count <- 0 # unused by predict(), formula only needs x on newdata
  spmod <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata_pois, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(0.7371, 0.0630))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.9334, 0.7299))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1088, 0.0824))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.5058, 0.5357))

  # autoregressive models

  spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car")
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1543, -0.1212))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2868))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1532, 0.1269))
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(1.0137))

  spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car", random = ~group, partition_factor = ~group)
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1481, -0.1226))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2807))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1442, 0.1275))
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(1.0147))

  # a different covariance family than the car cases above
  spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "sar")
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1337, -0.0751))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2369))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1362, 0.1282))
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(0.9492))

  spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly, spcov_type = "car")
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2769, 0.0317))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2277))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1119, 0.1049))
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(0.1932))

  spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly, spcov_type = "car", random = ~group, partition_factor = ~group)
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.2781, 0.0305))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(-0.2106))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1224, 0.1046))
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(0.1968))

  # a different GLM family than the Gamma cases above (a binary response)
  exdata_Mpoly_bin <- exdata_Mpoly
  exdata_Mpoly_bin$bin <- as.numeric(exdata_Mpoly_bin$y > median(exdata_Mpoly_bin$y, na.rm = TRUE))
  spmod <- spgautor(bin ~ x, exdata_Mpoly_bin, family = "binomial", spcov_type = "car")
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(0.0532, -0.1300))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), c(0.1424))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(11.4234, 0.2913))
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), c(0.7220))

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
  spmod <- splm(y ~ x, exdata,
    spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    anisotropy = TRUE, random = ~group, partition_factor = ~group, local = local_list
  )
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
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.2800, 0.3179))
  seed <- seed + n_skip

  set.seed(seed)
  spmod <- spglm(abs(y) ~ x, exdata,
    family = "Gamma", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    anisotropy = TRUE, random = ~group, partition_factor = ~group, local = local_list
  )
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE, local = local_list)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3349, 0.0437))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(-0.3043, -0.3981))
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.0938, 0.0852))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(0.2874, 0.2941))
})

test_that("exact output match for known covariance/dispersion parameters", {
  # characterization tests for the use_gloglik_known()/use_gloglik_known_anis()/
  # use_laploglik_known() dispatch paths
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

  randcov_initial_val <- randcov_initial(group = 1, known = "given")

  # use_gloglik_known() -- splm, isotropic
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, known = "given")
  spmod <- splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3482, -0.0718))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -135.6789)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.4177, 0.1149))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.5463, -0.8676))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(1.1912, 1.2612))

  spmod <- splm(y ~ x, exdata,
    spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord,
    random = ~group, randcov_initial = randcov_initial_val
  )
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3647, -0.0766))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -139.3669)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.6518, 0.1152))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.6278, -1.1097))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(1.2108, 1.2779))

  # use_gloglik_known_anis() -- splm, anisotropic
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, rotate = 2, scale = 0.5, known = "given")
  spmod <- splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3179, -0.0669))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -137.8189)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.3383, 0.1185))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.5665, -0.7861))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(1.2439, 1.2968))

  spmod <- splm(y ~ x, exdata,
    spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord,
    random = ~group, randcov_initial = randcov_initial_val
  )
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3312, -0.0693))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -141.4543)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.6039, 0.1187))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.6456, -1.0274))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(1.2652, 1.3168))

  # use_gloglik_known() -- spautor
  spcov_initial_val <- spcov_initial("car", de = 1, ie = 0, range = 0.5, extra = 1, known = "given")
  spmod <- spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val)
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1559, -0.1031))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -103.5713)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.0968, 0.0643))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), -0.1784)
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), 0.5133)

  spmod <- spautor(y ~ x, exdata_Mpoly,
    spcov_initial = spcov_initial_val,
    random = ~group, randcov_initial = randcov_initial_val
  )
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.1610, -0.1000))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -105.1860)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.7137, 0.0644))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), -0.2583)
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), 0.5178)

  # use_laploglik_known() -- spglm, isotropic (no prior coverage at all)
  exdata_pois <- exdata
  exdata_pois$count <- round(abs(exdata_pois$y) * 3)
  newexdata_pois <- newexdata
  newexdata_pois$count <- 0 # unused by predict(), formula only needs x on newdata

  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, known = "given")
  dispersion_initial_val <- dispersion_initial("poisson", dispersion = 1, known = "given")
  spmod <- spglm(count ~ x, exdata_pois,
    family = "poisson", spcov_initial = spcov_initial_val,
    dispersion_initial = dispersion_initial_val, xcoord = xcoord, ycoord = ycoord
  )
  preds <- predict(spmod, newdata = newexdata_pois, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(0.6479, 0.0949))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -298.7243)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.4312, 0.1392))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(0.9809, 0.8582))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(1.2078, 1.2783))

  spmod <- spglm(count ~ x, exdata_pois,
    family = "poisson", spcov_initial = spcov_initial_val,
    dispersion_initial = dispersion_initial_val, xcoord = xcoord, ycoord = ycoord,
    random = ~group, randcov_initial = randcov_initial_val
  )
  preds <- predict(spmod, newdata = newexdata_pois, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(0.6466, 0.0948))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -302.0582)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.6605, 0.1396))
  expect_equal(round(as.vector(preds$fit[1:2]), digits = 4), c(1.0206, 0.9238))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(1.2343, 1.3001))

  # use_laploglik_known() -- spgautor (no prior coverage at all)
  spcov_initial_val <- spcov_initial("car", de = 1, ie = 0, range = 0.5, extra = 1, known = "given")
  dispersion_initial_val <- dispersion_initial("Gamma", dispersion = 1, known = "given")
  spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly,
    spcov_initial = spcov_initial_val, dispersion_initial = dispersion_initial_val
  )
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3581, 0.0365))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -84.0287)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.1797, 0.1550))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), -0.3396)
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), 0.5843)

  spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly,
    spcov_initial = spcov_initial_val, dispersion_initial = dispersion_initial_val,
    random = ~group, randcov_initial = randcov_initial_val
  )
  preds <- predict(spmod, se.fit = TRUE)
  expect_equal(round(as.vector(coef(spmod)), digits = 4), c(-0.3546, 0.0327))
  expect_equal(round(as.numeric(logLik(spmod)), digits = 4), -85.3004)
  expect_equal(round(as.vector(sqrt(diag(vcov(spmod)))), digits = 4), c(0.7296, 0.1544))
  expect_equal(round(as.vector(preds$fit[1]), digits = 4), -0.2426)
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), 0.6059)
})

test_that("exact output match for prediction and confidence intervals", {
  # characterization tests for interval = "prediction"/"confidence"
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

  # splm, point-referenced Gaussian
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata, interval = "prediction")
  expect_equal(round(as.vector(preds[1:2, ]), digits = 4), c(0.6361, -0.9631, -0.8628, -2.5840, 2.1350, 0.6578))
  preds <- predict(spmod, newdata = newexdata, interval = "confidence")
  expect_equal(round(as.vector(preds[1:2, ]), digits = 4), c(-0.4454, -0.2885, -1.2085, -1.0614, 0.3177, 0.4844))

  # spautor, areal Gaussian
  spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car")
  preds <- predict(spmod, interval = "prediction")
  expect_equal(round(as.vector(preds[1, ]), digits = 4), c(-0.2868, -2.2737, 1.7000))
  preds <- predict(spmod, interval = "confidence")
  expect_equal(round(as.vector(preds[1, ]), digits = 4), c(-0.3428, -0.7981, 0.1125))

  # spglm, point-referenced GLM (poisson) -- interval = "prediction" with
  # se.fit and delta = TRUE together
  exdata_pois <- exdata
  exdata_pois$count <- round(abs(exdata_pois$y) * 3)
  newexdata_pois <- newexdata
  newexdata_pois$count <- 0 # unused by predict(), formula only needs x on newdata
  spmod <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata_pois, interval = "prediction", se.fit = TRUE, type = "response", delta = TRUE)
  expect_equal(round(as.vector(preds$fit[1:2, ]), digits = 4), c(2.5432, 2.0748, 0.9437, 0.7261, 6.8537, 5.9287))
  expect_equal(round(as.vector(preds$se.fit[1:2]), digits = 4), c(1.2864, 1.1115))
  preds <- predict(spmod, newdata = newexdata_pois, interval = "confidence", type = "response")
  expect_equal(round(as.vector(preds[1:2, ]), digits = 4), c(2.1841, 1.9076, 1.7255, 1.3776, 2.7646, 2.6415))

  # spgautor, areal GLM (Gamma)
  spmod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly, spcov_type = "car")
  preds <- predict(spmod, interval = "prediction", se.fit = TRUE, type = "response", delta = TRUE)
  expect_equal(round(as.vector(preds$fit[1, ]), digits = 4), c(0.7964, 0.5453, 1.1630))
  expect_equal(round(as.vector(preds$se.fit[1]), digits = 4), 0.1539)
  preds <- predict(spmod, interval = "confidence", type = "response")
  expect_equal(round(as.vector(preds[1, ]), digits = 4), c(0.7964, 0.5524, 1.1481))
})
