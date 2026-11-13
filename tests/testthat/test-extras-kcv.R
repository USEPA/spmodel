skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

test_that("local kcv reuses full-data coefficient covariance without theoretical refits", {
  n <- 24
  dat <- data.frame(
    y = 2 + (seq_len(n) %% 5), x = sin(seq_len(n)),
    cx = seq_len(n) / 6, cy = rep(0:3, 6), off = 0.05 * cos(seq_len(n))
  )
  groups <- rep(1:3, each = 8)
  folds <- rep(1:3, 8)
  params <- spcov_initial("exponential", de = 0.3, ie = 0.4, range = 3, known = "given")
  fit_args <- list(
    formula = y ~ x + offset(off), data = dat, spcov_initial = params,
    xcoord = "cx", ycoord = "cy", estmethod = "ml",
    local = list(index = groups, var_adjust = "theoretical")
  )
  models <- list(
    splm = do.call(splm, c(fit_args, list(ddf = "asymptotic"))),
    spglm = do.call(spglm, c(fit_args, list(family = "poisson")))
  )
  local_mocked_bindings(
    get_W_ij = function(...) stop("Unexpected theoretical variance adjustment"),
    .package = "spmodel"
  )

  for (model_type in names(models)) {
    object <- models[[model_type]]
    original_vcov <- object$vcov
    original_beta <- coef(object)
    pred_local <- list(method = "distance", size = 6, parallel = FALSE)
    expected_fit <- expected_se <- unadjusted_se <- numeric(n)
    beta_changed <- FALSE
    for (fold in unique(folds)) {
      rows <- which(folds == fold)
      args <- fit_args
      args$data$y[rows] <- NA
      args$local <- list(index = groups[-rows], var_adjust = "none")
      if (model_type == "splm") {
        args$ddf <- "asymptotic"
      } else {
        args$family <- "poisson"
      }
      refit <- do.call(model_type, args)
      beta_changed <- beta_changed || max(abs(coef(refit) - original_beta)) > 1e-6
      unadjusted_se[rows] <- predict(refit, se.fit = TRUE, local = pred_local)$se.fit
      if (model_type == "splm") {
        refit$vcov$fixed <- vcov(object)
      } else {
        refit$vcov$fixed <- list(
          corrected = vcov(object), uncorrected = vcov(object, var_correct = FALSE)
        )
      }
      reference <- predict(refit, se.fit = TRUE, local = pred_local)
      expected_fit[rows] <- reference$fit
      expected_se[rows] <- reference$se.fit
    }
    expect_true(beta_changed)
    expect_gt(max(abs(expected_se - unadjusted_se)), 1e-6)
    fit_local <- get_kcv_estimation_local(object, which(folds == 1), get_local_list_prediction(pred_local))
    expect_identical(fit_local$var_adjust, "none")
    expect_equal(fit_local$index, groups[folds != 1])
    actual <- kcv(object, folds_index = folds, local = pred_local, cv_predict = TRUE, se.fit = TRUE)
    expect_equal(actual$cv_predict, expected_fit, tolerance = 1e-8)
    expect_equal(actual$se.fit, expected_se, tolerance = 1e-8)
    without_se <- kcv(object, folds_index = folds, local = pred_local, cv_predict = TRUE)
    expect_equal(without_se$cv_predict, expected_fit, tolerance = 1e-8)
    expect_identical(object$vcov, original_vcov)
    expect_identical(coef(object), original_beta)
    if (model_type == "splm") {
      intervals <- kcv(object, folds_index = folds, local = pred_local, interval = "prediction", level = 0.8)
      coverage <- mean(abs(dat$y - expected_fit) <= qnorm(0.9) * expected_se)
      expect_equal(intervals$cover.8, coverage)
    }
  }
})

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

  # Prediction-neighborhood methods are not valid fitting methods; kcv()
  # must translate them before refitting each fold.
  folds <- rep(1:5, length.out = n)
  kc_covariance <- kcv(spmod1,
    folds_index = folds, cv_predict = TRUE, se.fit = TRUE,
    local = list(method = "covariance", size = 12, parallel = FALSE)
  )
  expect_true(all(is.finite(kc_covariance$cv_predict)))
  expect_true(all(is.finite(kc_covariance$se.fit)))

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

  folds <- rep(1:5, length.out = n)
  kc_distance <- kcv(spmod1,
    folds_index = folds, cv_predict = TRUE, se.fit = TRUE,
    local = list(method = "distance", size = 12, parallel = FALSE)
  )
  expect_true(all(is.finite(kc_distance$cv_predict)))
  expect_true(all(is.finite(kc_distance$se.fit)))

  prediction_local <- get_local_list_prediction(list(method = "covariance", size = 12, parallel = FALSE))
  fit_local <- get_kcv_estimation_local(spmod1, fold_rows = 1:4, prediction_local)
  expect_equal(fit_local$method, "kmeans")
  expect_equal(fit_local$size, 12)
  spmod1$local_index <- rep(1:3, length.out = n)
  fit_local_indexed <- get_kcv_estimation_local(spmod1, fold_rows = 1:4, prediction_local)
  expect_equal(fit_local_indexed$index, spmod1$local_index[-(1:4)])

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
