# SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
load(system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

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

  smod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "penta", estmethod = "reml")
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

