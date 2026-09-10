skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

# SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

test_that("kcv works splm point data", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  spmod1 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  n <- spmod1$n

  # default k (5), stats-only return
  expect_s3_class(kcv(spmod1), "tbl_df")
  expect_named(kcv(spmod1), c("bias", "MSPE", "RMSPE", "cor2"))

  # cv_predict/se.fit return shape and length
  kc <- kcv(spmod1, cv_predict = TRUE, se.fit = TRUE)
  expect_type(kc, "list")
  expect_named(kc, c("stats", "cv_predict", "se.fit"))
  expect_length(kc$cv_predict, n)
  expect_length(kc$se.fit, n)
  expect_false(anyNA(kc$cv_predict))
  expect_false(anyNA(kc$se.fit))

  # k = n delegates to (and exactly matches) loocv()
  lo <- loocv(spmod1, cv_predict = TRUE, se.fit = TRUE)
  kn <- kcv(spmod1, k = n, cv_predict = TRUE, se.fit = TRUE)
  expect_equal(lo$stats, kn$stats)
  expect_equal(lo$cv_predict, kn$cv_predict)
  expect_equal(lo$se.fit, kn$se.fit)

  # local (big data) path still runs and returns sensible output even on
  # small data, when explicitly requested
  kc_local <- kcv(spmod1, local = TRUE, cv_predict = TRUE, se.fit = TRUE)
  expect_type(kc_local, "list")
  expect_length(kc_local$cv_predict, n)
  expect_false(anyNA(kc_local$cv_predict))
  expect_false(anyNA(kc_local$se.fit))

  # interval = "prediction": adds cover.95 to $stats, coverage matches
  # an independent manual computation from cv_predict/se.fit, and interval =
  # "none" (the default) is unaffected
  lo_pred <- loocv(spmod1, interval = "prediction")
  expect_named(lo_pred, c("bias", "MSPE", "RMSPE", "cor2", "cover.95"))
  expect_true(lo_pred$cover.95 >= 0 && lo_pred$cover.95 <= 1)

  y <- model.response(model.frame(spmod1))
  tstar <- qnorm(0.975)
  lwr <- lo$cv_predict - tstar * lo$se.fit
  upr <- lo$cv_predict + tstar * lo$se.fit
  expect_equal(lo_pred$cover.95, mean(y >= lwr & y <= upr))

  kc_pred <- kcv(spmod1, k = n, interval = "prediction")
  expect_equal(kc_pred, lo_pred)

  lo_default <- loocv(spmod1)
  expect_named(lo_default, c("bias", "MSPE", "RMSPE", "cor2"))

  expect_error(loocv(spmod1, interval = "confidence"))
  expect_error(kcv(spmod1, interval = "confidence"))

  # fold sizes are (approximately) equal
  set.seed(1)
  fold_id <- get_kcv_folds(5, n)
  expect_equal(length(fold_id), n)
  expect_lte(diff(range(table(fold_id))), 1)

  # errors
  expect_error(kcv(spmod1, k = n + 1), "k cannot exceed")
  expect_error(kcv(spmod1, k = 1), "k must be a single whole number")
  expect_error(kcv(spmod1, k = 2.5), "k must be a single whole number")
})

test_that("kcv works spautor polygon data", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  spmod1 <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")
  n <- spmod1$n

  expect_s3_class(kcv(spmod1), "tbl_df")
  expect_named(kcv(spmod1), c("bias", "MSPE", "RMSPE", "cor2"))

  kc <- kcv(spmod1, cv_predict = TRUE, se.fit = TRUE)
  expect_type(kc, "list")
  expect_length(kc$cv_predict, n)
  expect_length(kc$se.fit, n)
  expect_false(anyNA(kc$cv_predict))
  expect_false(anyNA(kc$se.fit))

  lo <- loocv(spmod1, cv_predict = TRUE, se.fit = TRUE)
  kn <- kcv(spmod1, k = n, cv_predict = TRUE, se.fit = TRUE)
  expect_equal(lo$stats, kn$stats)
  expect_equal(lo$cv_predict, kn$cv_predict)
  expect_equal(lo$se.fit, kn$se.fit)

  # interval = "prediction"
  lo_pred <- loocv(spmod1, interval = "prediction")
  expect_named(lo_pred, c("bias", "MSPE", "RMSPE", "cor2", "cover.95"))
  expect_true(lo_pred$cover.95 >= 0 && lo_pred$cover.95 <= 1)

  y <- model.response(model.frame(spmod1))
  tstar <- qnorm(0.975)
  lwr <- lo$cv_predict - tstar * lo$se.fit
  upr <- lo$cv_predict + tstar * lo$se.fit
  expect_equal(lo_pred$cover.95, mean(y >= lwr & y <= upr))

  kc_pred <- kcv(spmod1, k = n, interval = "prediction")
  expect_equal(kc_pred, lo_pred)

  expect_error(kcv(spmod1, k = n + 1), "k cannot exceed")
  expect_error(kcv(spmod1, k = 1), "k must be a single whole number")
})

test_that("kcv works spglm point data", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  spmod1 <- spglm(abs(y) ~ x, "Gamma", exdata, spcov_type = "exponential", xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml")
  n <- spmod1$n

  expect_s3_class(kcv(spmod1), "tbl_df")
  expect_named(kcv(spmod1), c("bias", "MSPE", "RMSPE"))

  kc <- kcv(spmod1, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE)
  expect_type(kc, "list")
  expect_length(kc$cv_predict, n)
  expect_length(kc$se.fit, n)
  expect_false(anyNA(kc$cv_predict))
  expect_false(anyNA(kc$se.fit))

  lo <- loocv(spmod1, cv_predict = TRUE, se.fit = TRUE)
  kn <- kcv(spmod1, k = n, cv_predict = TRUE, se.fit = TRUE)
  expect_equal(lo$stats, kn$stats)
  expect_equal(lo$cv_predict, kn$cv_predict)
  expect_equal(lo$se.fit, kn$se.fit)

  kc_local <- kcv(spmod1, local = TRUE, cv_predict = TRUE)
  expect_false(anyNA(kc_local$cv_predict))

  expect_error(kcv(spmod1, k = n + 1), "k cannot exceed")
  expect_error(kcv(spmod1, k = 1), "k must be a single whole number")
})

test_that("kcv works spgautor polygon data", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  spmod1 <- spgautor(abs(y) ~ x, family = Gamma, exdata_poly, spcov_type = "car", estmethod = "reml")
  n <- spmod1$n

  expect_s3_class(kcv(spmod1), "tbl_df")
  expect_named(kcv(spmod1), c("bias", "MSPE", "RMSPE"))

  kc <- kcv(spmod1, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE)
  expect_type(kc, "list")
  expect_length(kc$cv_predict, n)
  expect_length(kc$se.fit, n)
  expect_false(anyNA(kc$cv_predict))
  expect_false(anyNA(kc$se.fit))

  lo <- loocv(spmod1, cv_predict = TRUE, se.fit = TRUE)
  kn <- kcv(spmod1, k = n, cv_predict = TRUE, se.fit = TRUE)
  expect_equal(lo$stats, kn$stats)
  expect_equal(lo$cv_predict, kn$cv_predict)
  expect_equal(lo$se.fit, kn$se.fit)

  expect_error(kcv(spmod1, k = n + 1), "k cannot exceed")
  expect_error(kcv(spmod1, k = 1), "k must be a single whole number")
})
