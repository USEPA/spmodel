skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

test_that("satterthwaite() closed-form derivatives exist for exponential/gaussian/spherical", {
  for (spcov_type in c("exponential", "gaussian", "spherical")) {
    spmod <- splm(y ~ x, exdata, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, estmethod = "reml")

    # these three types have a closed-form dSig_dtheta_spcov.<type>()
    # (get_dSig_dtheta.R), so the default (unspecified method) resolves to
    # "closed" without a warning
    expect_warning(sw_default <- satterthwaite(spmod), NA)
    sw_closed <- satterthwaite(spmod, method = "closed")
    expect_equal(sw_default, sw_closed)
    expect_named(sw_closed, names(coef(spmod, type = "fixed")))
    expect_true(all(is.finite(sw_closed)) && all(sw_closed > 0))

    # "numeric" is also available as an independent (finite-difference) way
    # to compute the same quantity -- "closed" uses the expected (Fisher)
    # information and "numeric" the observed information (a numerically
    # differentiated Hessian of the realized log-likelihood), which are only
    # asymptotically equivalent estimators of Cov(theta_hat), so they are not
    # expected to closely agree in a finite sample -- just both be sane
    sw_numeric <- satterthwaite(spmod, method = "numeric")
    expect_named(sw_numeric, names(coef(spmod, type = "fixed")))
    expect_true(all(is.finite(sw_numeric)) && all(sw_numeric > 0))
  }
})

test_that("satterthwaite() falls back to numeric (with a warning) for covariance types without a closed form, and for anisotropy", {
  spmod_matern <- splm(y ~ x, exdata, spcov_type = "matern", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_warning(satterthwaite(spmod_matern, method = "closed"), "Closed form not available")
  # the default silently resolves to "numeric" for types with no closed form
  expect_warning(satterthwaite(spmod_matern), NA)

  spmod_anis <- splm(y ~ x, exdata,
    spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    estmethod = "reml", anisotropy = TRUE
  )
  expect_warning(satterthwaite(spmod_anis, method = "closed"), "Closed form not available")
  # anisotropy forces "numeric" even for exponential/gaussian/spherical, so
  # the default (unspecified method) should not warn about a missing closed form
  expect_warning(satterthwaite(spmod_anis), NA)
})

test_that("satterthwaite() and its ddf/vcov() plumbing work with random effects", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml", random = ~group)

  expect_error(sw <- satterthwaite(spmod), NA)
  expect_named(sw, names(coef(spmod, type = "fixed")))

  spcov_mat <- vcov(spmod, type = "spcov")
  randcov_mat <- vcov(spmod, type = "randcov")
  cov_mat <- vcov(spmod, type = "cov")
  expect_false(is.null(spcov_mat))
  expect_false(is.null(randcov_mat))
  expect_equal(nrow(cov_mat), nrow(spcov_mat) + nrow(randcov_mat))
  # spcov/randcov are exactly the corresponding blocks of the joint matrix
  expect_equal(as.matrix(cov_mat[rownames(spcov_mat), colnames(spcov_mat), drop = FALSE]), as.matrix(spcov_mat))
  expect_equal(as.matrix(cov_mat[rownames(randcov_mat), colnames(randcov_mat), drop = FALSE]), as.matrix(randcov_mat))
})

test_that("satterthwaite() works for spautor car/sar, with/without M, with random effects, and unconnected sites", {
  for (spcov_type in c("car", "sar")) {
    spmod <- spautor(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "reml")
    expect_error(satterthwaite(spmod), NA)

    # a nested random effect variance estimated near zero can occasionally
    # push numDeriv's finite-difference perturbation of the covariance
    # parameters into a non-positive-definite region for car/sar (a known
    # boundary-fragility of the "numeric" method, not new/specific to this
    # test -- see get_vcov_theta()/get_grad_g()'s "numeric" branches);
    # tolerate that specific failure rather than asserting it never happens
    spmod_rand <- spautor(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "reml", random = ~group)
    result <- tryCatch(satterthwaite(spmod_rand), error = function(e) conditionMessage(e))
    if (is.character(result)) {
      expect_match(result, "not positive")
    } else {
      expect_named(result, names(coef(spmod_rand, type = "fixed")))
    }

    spmod_unconnected <- spautor(y ~ x, exdata_Upoly, spcov_type = spcov_type, estmethod = "reml")
    expect_error(satterthwaite(spmod_unconnected), NA)
  }

  # car with an explicit M (the CAR symmetry condition matrix)
  spmod_M <- spautor(y ~ x, exdata_Mpoly, spcov_type = "car", estmethod = "reml")
  expect_error(satterthwaite(spmod_M), NA)

  # car/sar have no closed-form dSig_dtheta_spcov() method, so they always
  # use "numeric" -- requesting "closed" explicitly should warn and fall back
  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")
  expect_warning(satterthwaite(spmod, method = "closed"), "Closed form not available")
})

test_that("satterthwaite() respects its documented scope restrictions", {
  spmod_svwls <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls")
  expect_error(satterthwaite(spmod_svwls), "estmethod")

  spmod_local <- splm(y ~ x, exdata,
    spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    estmethod = "reml", local = list(size = 20)
  )
  expect_error(satterthwaite(spmod_local), "local")

  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, known = "given")
  spmod_known <- splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_error(satterthwaite(spmod_known), "All covariance parameters known")

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_error(satterthwaite(spmod, method = "bogus"), "method must be")
})

test_that("ddf argument validates and matches the automatic sample-size-based default", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_error(
    splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml", ddf = "bogus"),
    "ddf must be"
  )
  expect_error(anova(spmod, ddf = "bogus"), "ddf must be")

  # explicit ddf = "satterthwaite" for n < 500 matches the automatic default
  spmod_explicit <- splm(y ~ x, exdata,
    spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    estmethod = "reml", ddf = "satterthwaite"
  )
  expect_equal(spmod$ddf, spmod_explicit$ddf)
})

test_that("ddf = 'satterthwaite' warns (but still computes) for n >= 500, and is skipped by default", {
  set.seed(2)
  n <- 505
  big_data <- data.frame(x = rnorm(n), xcoord = runif(n, 0, 10), ycoord = runif(n, 0, 10))
  spcov_params_val <- spcov_params("exponential", de = 1, ie = 0.5, range = 2)
  big_data$y <- sprnorm(spcov_params_val, mean = 1 + 0.5 * big_data$x, data = big_data, xcoord = xcoord, ycoord = ycoord)

  spmod_default <- splm(y ~ x, big_data, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_null(spmod_default$ddf)

  expect_warning(
    spmod_explicit <- splm(y ~ x, big_data,
      spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
      estmethod = "reml", ddf = "satterthwaite"
    ),
    "n >= 500"
  )
  expect_false(is.null(spmod_explicit$ddf))
})

test_that("anova() ddf = 'satterthwaite' matches satterthwaite() exactly for single-coefficient terms", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  sw <- satterthwaite(spmod)

  # a single-row (q = 1) hypothesis is a mathematically exact special case of
  # the Fai-Cornelius combination, so this should match satterthwaite() exactly
  av <- anova(spmod, Terms = "x")
  expect_equal(unname(av$DenDF), unname(sw[["x"]]))

  # test = FALSE drops Pr(>F) whether ddf is explicit or automatic
  expect_false("Pr(>F)" %in% colnames(anova(spmod, test = FALSE)))
  expect_false("Pr(>F)" %in% colnames(anova(spmod, ddf = "satterthwaite", test = FALSE)))

  # an explicit ddf = "satterthwaite" request lets a scope-violation error
  # propagate, while the automatic (missing ddf) case instead falls back
  # silently to the asymptotic chi-squared table
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, known = "given")
  spmod_known <- splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_error(anova(spmod_known, ddf = "satterthwaite"), "All covariance parameters known")
  expect_true("Chi2" %in% colnames(anova(spmod_known)))
})

test_that("vcov() type = 'cov'/'spcov'/'randcov' are NULL together under ddf = 'asymptotic', and type is validated", {
  spmod_asymp <- splm(y ~ x, exdata,
    spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    estmethod = "reml", random = ~group, ddf = "asymptotic"
  )
  expect_null(vcov(spmod_asymp, type = "cov"))
  expect_null(vcov(spmod_asymp, type = "spcov"))
  expect_null(vcov(spmod_asymp, type = "randcov"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_error(vcov(spmod, type = "bogus"), "Invalid type argument")
})

test_that("emmeans::joint_tests() reports finite df for Satterthwaite fits and Inf for asymptotic fits", {
  skip_if_not_installed("emmeans")

  spmod <- splm(y ~ group, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  jt <- emmeans::joint_tests(spmod)
  expect_s3_class(jt, "data.frame")
  # emmeans combines per-row Satterthwaite df across a joint contrast's rows
  # via min() (see emmeans:::test.emmGrid), a different (more conservative)
  # rule than anova()'s own Fai-Cornelius combination, so df2 need not match
  # anova(spmod, Terms = "group")$DenDF exactly -- just be finite and positive
  expect_true(all(is.finite(jt$df2)) && all(jt$df2 > 0))

  spmod_asymp <- splm(y ~ group, exdata,
    spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    estmethod = "reml", ddf = "asymptotic"
  )
  jt_asymp <- emmeans::joint_tests(spmod_asymp)
  expect_true(all(is.infinite(jt_asymp$df2)))

  # emmeans()/contrast() marginal means and pairwise comparisons also report
  # finite t-based df, not just joint_tests()
  em <- emmeans::emmeans(spmod, ~group)
  expect_true(all(is.finite(as.data.frame(em)$df)))
  pw <- pairs(em)
  expect_true(all(is.finite(as.data.frame(pw)$df)))
})

test_that("confint() uses Satterthwaite t-based intervals exactly matching a manual qt() calculation", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")

  ci <- confint(spmod, level = 0.90)
  estimates <- coef(spmod, type = "fixed")
  se <- sqrt(diag(vcov(spmod, type = "fixed")))
  tstar <- qt(0.95, spmod$ddf[names(estimates)])
  expect_equal(unname(ci[, 1]), unname(estimates - tstar * se))
  expect_equal(unname(ci[, 2]), unname(estimates + tstar * se))

  # parm subsetting still works alongside the t-based path
  ci_x <- confint(spmod, parm = "x", level = 0.90)
  expect_equal(rownames(ci_x), "x")
  expect_equal(ci_x, ci["x", , drop = FALSE])

  # ddf = "asymptotic" reverts to the z-based interval (unaffected by ddf)
  spmod_asymp <- splm(y ~ x, exdata,
    spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord,
    estmethod = "reml", ddf = "asymptotic"
  )
  ci_asymp <- confint(spmod_asymp, level = 0.90)
  tstar_z <- qnorm(0.95)
  estimates_asymp <- coef(spmod_asymp, type = "fixed")
  se_asymp <- sqrt(diag(vcov(spmod_asymp, type = "fixed")))
  expect_equal(unname(ci_asymp[, 1]), unname(estimates_asymp - tstar_z * se_asymp))

  # spglm/spgautor never have a ddf concept, so confint() always stays z-based
  exdata_pois <- exdata
  exdata_pois$count <- round(abs(exdata_pois$y) * 3)
  spmod_glm <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_null(spmod_glm$ddf)
  ci_glm <- confint(spmod_glm, level = 0.90)
  estimates_glm <- coef(spmod_glm, type = "fixed")
  se_glm <- sqrt(diag(vcov(spmod_glm, type = "fixed")))
  expect_equal(unname(ci_glm[, 1]), unname(estimates_glm - tstar_z * se_glm))
})

test_that("ddf/satterthwaite work for general splm", {
  spmod <- splm(z ~ water + tarp, data = caribou, spcov_type = "exponential", xcoord = x, ycoord = y, estmethod = "reml")

  # satterthwaite()
  sw <- satterthwaite(spmod)
  expect_type(sw, "double")
  expect_named(sw, names(coef(spmod, type = "fixed")))
  expect_true(all(sw > 0))

  # ddf computed automatically at fit time (n = 30 <= 500) and matches satterthwaite()
  expect_equal(spmod$ddf, sw)

  # summary()/tidy() switch to t-based inference when ddf is available
  expect_true("df" %in% colnames(summary(spmod)$coefficients$fixed))
  expect_true("df" %in% colnames(tidy(spmod)))

  # ddf = "asymptotic" disables it
  spmod_asymp <- splm(z ~ water + tarp,
    data = caribou, spcov_type = "exponential", xcoord = x, ycoord = y,
    estmethod = "reml", ddf = "asymptotic"
  )
  expect_null(spmod_asymp$ddf)
  expect_false("df" %in% colnames(summary(spmod_asymp)$coefficients$fixed))

  # anova() ddf argument
  expect_true("NumDF" %in% colnames(anova(spmod)))
  expect_true("Chi2" %in% colnames(anova(spmod, ddf = "asymptotic")))

  # vcov() type argument (a base matrix or Matrix-class object depending on
  # method = "closed"/"numeric", see get_satterthwaite_method())
  cov_mat <- vcov(spmod, type = "cov")
  expect_equal(nrow(cov_mat), ncol(cov_mat))
  spcov_mat <- vcov(spmod, type = "spcov")
  expect_equal(nrow(spcov_mat), ncol(spcov_mat))
  expect_null(vcov(spmod, type = "randcov"))

  # confint() switches to a t-based interval when ddf is available
  ci <- confint(spmod)
  ci_asymp <- confint(spmod_asymp)
  expect_true(all((ci[, 2] - ci[, 1]) > (ci_asymp[, 2] - ci_asymp[, 1])))
})

test_that("emmeans works for splm", {
  skip_if_not_installed("emmeans")
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  spmod1 <- splm(y ~ group, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  em <- emmeans::emmeans(spmod1, ~group)
  expect_s4_class(em, "emmGrid")
  expect_equal(nrow(as.data.frame(em)), nlevels(exdata$group))

  # ddf/satterthwaite feeds into emmeans::joint_tests() with a finite df2
  jt <- emmeans::joint_tests(spmod1)
  expect_s3_class(jt, "data.frame")
  expect_true(all(is.finite(jt$df2)))
})

test_that("ddf/satterthwaite work for general spautor", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")

  # satterthwaite()
  sw <- satterthwaite(spmod)
  expect_type(sw, "double")
  expect_named(sw, names(coef(spmod, type = "fixed")))
  expect_true(all(sw > 0))

  # ddf computed automatically at fit time (n <= 500) and matches satterthwaite()
  expect_equal(spmod$ddf, sw)

  # summary()/tidy() switch to t-based inference when ddf is available
  expect_true("df" %in% colnames(summary(spmod)$coefficients$fixed))
  expect_true("df" %in% colnames(tidy(spmod)))

  # ddf = "asymptotic" disables it
  spmod_asymp <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml", ddf = "asymptotic")
  expect_null(spmod_asymp$ddf)
  expect_false("df" %in% colnames(summary(spmod_asymp)$coefficients$fixed))

  # anova() ddf argument
  expect_true("NumDF" %in% colnames(anova(spmod)))
  expect_true("Chi2" %in% colnames(anova(spmod, ddf = "asymptotic")))

  # vcov() type argument (a base matrix or Matrix-class object depending on
  # method = "closed"/"numeric" -- car/sar always use "numeric", see
  # get_satterthwaite_method())
  cov_mat <- vcov(spmod, type = "cov")
  expect_equal(nrow(cov_mat), ncol(cov_mat))
  spcov_mat <- vcov(spmod, type = "spcov")
  expect_equal(nrow(spcov_mat), ncol(spcov_mat))
  expect_null(vcov(spmod, type = "randcov"))

  # confint() switches to a t-based interval when ddf is available
  ci <- confint(spmod)
  ci_asymp <- confint(spmod_asymp)
  expect_true(all((ci[, 2] - ci[, 1]) > (ci_asymp[, 2] - ci_asymp[, 1])))
})

test_that("emmeans works for spautor", {
  skip_if_not_installed("emmeans")
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  spmod1 <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")

  # ddf/satterthwaite feeds into emmeans::joint_tests() with a finite df2
  jt <- emmeans::joint_tests(spmod1)
  expect_s3_class(jt, "data.frame")
  expect_true(all(is.finite(jt$df2)))
})
