skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

exdata_pois <- exdata
exdata_pois$count <- round(abs(exdata_pois$y) * 3)

# conditional() is only implemented for splm()/spglm(): for spautor()/spgautor(),
# the conditional simulation locations would need to be known ahead of fitting
# (they affect the areal neighborhood structure), so a base-and-block big-data
# approximation like the one used here is not feasible
# there is intentionally
# no conditional.spautor()/conditional.spgautor() method.
test_that("conditional() has no method for spautor/spgautor (out of scope by design)", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  spmod_autor <- spautor(y ~ x, exdata_poly, spcov_type = "car")
  expect_error(conditional(spmod_autor), "no applicable method")
})

test_that("conditional() validates output/type and the newdata/object$newdata requirement", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_error(conditional(spmod, newdata = newexdata, output = "bogus"), "output must be")

  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_error(conditional(spmod_g, newdata = newexdata, output = "bogus"), "output must be")
  expect_error(conditional(spmod_g, newdata = newexdata, type = "bogus"), "should be one of")

  # no newdata argument and no missing observations to fall back on
  expect_error(conditional(spmod), "No missing data to predict")
})

test_that("output argument returns the documented shapes for splm() and spglm()", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  samples <- 20

  cond_newdata <- conditional(spmod, newdata = newexdata, output = "newdata", samples = samples)
  expect_equal(dim(cond_newdata), c(NROW(newexdata), samples))

  cond_beta <- conditional(spmod, newdata = newexdata, output = "beta", samples = samples)
  expect_equal(dim(cond_beta), c(length(coef(spmod)), samples))

  cond_object <- conditional(spmod, newdata = newexdata, output = "object", samples = samples)
  expect_equal(dim(cond_object), c(spmod$n, samples))
  # output = "object" just replicates the observed response, so every column
  # (simulation draw) is identical and equal to the observed y
  expect_true(all(apply(cond_object, 1, function(row) length(unique(row)) == 1)))
  expect_equal(cond_object[, 1], unname(model.response(model.frame(spmod))))

  cond_multi <- conditional(spmod, newdata = newexdata, output = c("newdata", "beta"), samples = samples)
  expect_named(cond_multi, c("newdata", "beta"))

  cond_all <- conditional(spmod, newdata = newexdata, output = "all", samples = samples)
  expect_named(cond_all, c("newdata", "beta", "object"))

  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  cond_object_g <- conditional(spmod_g, newdata = newexdata, output = "object", samples = samples)
  expect_equal(dim(cond_object_g), c(spmod_g$n, samples))
  # for spglm(), "object" replicates the fitted link-scale latent process w,
  # not the observed response y
  expect_equal(cond_object_g[, 1], unname(fitted(spmod_g, type = "link")))
})

test_that("type argument controls the returned scale for spglm()", {
  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  samples <- 30

  # the type-scale conversion is the last step applied to an otherwise
  # identical sequence of random draws, so matching seeds makes the
  # link/response comparison below an exact check, not just approximate
  set.seed(101)
  cond_link <- conditional(spmod_g, newdata = newexdata, type = "link", samples = samples)
  set.seed(101)
  cond_response <- conditional(spmod_g, newdata = newexdata, type = "response", samples = samples)
  cond_new <- conditional(spmod_g, newdata = newexdata, type = "new", samples = samples)

  expect_true(all(cond_response >= 0))
  expect_true(all(cond_new >= 0))
  # a "new" draw for a Poisson response is a simulated count
  expect_equal(cond_new, round(cond_new))
  # the response-scale mean is exp() of the link-scale draw here (log link,
  # no dispersion adjustment for a Poisson response)
  expect_equal(cond_response, exp(cond_link))
})

test_that("conditional.spglm() newdata_size controls the size of simulated binomial draws", {
  exdata_bin <- exdata
  exdata_bin$bern <- rbinom(NROW(exdata_bin), size = 1, prob = 0.5)
  spmod_bin <- spglm(bern ~ x, family = binomial, data = exdata_bin, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential")

  cond_default <- conditional(spmod_bin, newdata = newexdata, type = "new", samples = 30)
  expect_true(all(cond_default %in% c(0, 1)))

  cond_size5 <- conditional(spmod_bin, newdata = newexdata, type = "new", samples = 30, newdata_size = rep(5, NROW(newexdata)))
  expect_true(all(cond_size5 >= 0 & cond_size5 <= 5))
  expect_true(any(cond_size5 > 1)) # exercising sizes above the default of 1
})

test_that("conditional() respects offset() identically to an equivalent shifted-response model", {
  exdata_off <- exdata
  exdata_off$offset <- 2
  newexdata_off <- newexdata
  newexdata_off$offset <- 2

  spmod_off <- splm(y ~ x + offset(offset), exdata_off, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  spmod_no_off <- splm(I(y - 2) ~ x, exdata_off, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  set.seed(3)
  cond_off <- conditional(spmod_off, newdata = newexdata_off, samples = 50)
  set.seed(3)
  cond_no_off <- conditional(spmod_no_off, newdata = newexdata_off, samples = 50)
  # matched seeds + an additive offset should reproduce the same draws shifted
  # by exactly the offset (2), since the offset is added back on deterministically
  expect_equal(cond_off, cond_no_off + 2)
})

test_that("conditional() works with random effects, anisotropy, a partition factor, and a pure-nugget covariance", {
  samples <- 20

  spmod_rand <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~group)
  expect_error(cond_rand <- conditional(spmod_rand, newdata = newexdata, samples = samples), NA)
  expect_true(all(is.finite(cond_rand)))

  spmod_anis <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, anisotropy = TRUE)
  expect_error(cond_anis <- conditional(spmod_anis, newdata = newexdata, samples = samples), NA)
  expect_true(all(is.finite(cond_anis)))

  spmod_pf <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, partition_factor = ~group)
  expect_error(cond_pf <- conditional(spmod_pf, newdata = newexdata, samples = samples), NA)
  expect_true(all(is.finite(cond_pf)))

  # spcov_type = "none": a pure-nugget model, exercised because
  # get_conditional_new_from_base_adjust() takes a diagonal-matrix shortcut
  # (element-wise sqrt() instead of chol()) in this case
  spmod_nugget <- splm(y ~ x, exdata, spcov_type = "none", xcoord = xcoord, ycoord = ycoord)
  expect_error(cond_nugget <- conditional(spmod_nugget, newdata = newexdata, samples = samples), NA)
  expect_true(all(is.finite(cond_nugget)))
})

test_that("conditional() SD is close to predict() se.fit under ordinary (non-extreme) conditions", {
  # a broad sanity check, not a tight regression test (see the dedicated
  # regression tests below for that): conditional() draws should have
  # roughly the same spread as predict()'s analytic standard error
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  set.seed(11)
  cond <- conditional(spmod, newdata = newexdata, samples = 3000)
  ratio <- apply(cond, 1, sd) / preds$se.fit
  expect_true(all(ratio > 0.8 & ratio < 1.25))
  # simulated means should track predict()'s fitted values (no bias); an
  # absolute threshold is used since the fitted values straddle zero, which
  # would make a relative tolerance behave inconsistently across rows
  expect_true(all(abs(apply(cond, 1, mean) - preds$fit) < 0.1))

  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  preds_g <- predict(spmod_g, newdata = newexdata, type = "link", se.fit = TRUE)
  set.seed(12)
  cond_g <- conditional(spmod_g, newdata = newexdata, samples = 3000)
  ratio_g <- apply(cond_g, 1, sd) / preds_g$se.fit
  expect_true(all(ratio_g > 0.8 & ratio_g < 1.25))
  expect_true(all(abs(apply(cond_g, 1, mean) - preds_g$fit) < 0.1))
})

test_that("conditional.splm() no longer double-counts fixed effect uncertainty (regression test)", {
  # get_conditional_new_from_base_adjust() used to add an analytic
  # H %*% cov_betahat %*% t(H) term on top of conditional draws that were
  # already built from simulated beta draws and double-counting fixed effect
  # uncertainty. That term scales with the newdata design point x0 (via H),
  # so a covariate value far outside the observed range makes it dominate the
  # (otherwise small, since this point is spatially coincident with an
  # observed location) kriging variance and turning a subtle miscalibration
  # into an unmistakable one if the bug ever returns.
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_true(50 > max(exdata$x) + 10) # confirm x = 50 is well outside the observed range

  newdata_extreme <- newexdata[1, , drop = FALSE]
  newdata_extreme$x <- 50
  newdata_extreme$xcoord <- exdata$xcoord[1]
  newdata_extreme$ycoord <- exdata$ycoord[1]

  preds <- predict(spmod, newdata = newdata_extreme, se.fit = TRUE)
  set.seed(21)
  cond <- conditional(spmod, newdata = newdata_extreme, samples = 3000)
  ratio <- sd(cond) / preds$se.fit
  # the (fixed) implementation should land close to 1; the removed
  # double-counting term would have pushed this well above 1.15 here
  expect_true(ratio > 0.85 && ratio < 1.15)
})

test_that("conditional.spglm() draws are unchanged from a known-good seeded snapshot (regression test)", {
  # a deterministic canary for conditional.spglm()'s random-draw structure:
  # this pins down the exact sequence of rnorm() draws consumed. Both fixed
  # bugs (the removed w-resimulation block and the removed analytic
  # H %*% cov_betahat %*% t(H) term) changed how many/which random draws were
  # consumed, so reintroducing either one would change every value below.
  set.seed(42)
  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  set.seed(42)
  cond_g <- conditional(spmod_g, newdata = newexdata[1:2, ], samples = 4)

  expect_equal(
    round(as.vector(cond_g), 4),
    c(1.9616, 0.8677, 1.6446, 1.8974, 0.2404, 0.6272, 0.9122, 1.1972)
  )
})

test_that("conditional.splm() draws are unchanged from a known-good seeded snapshot (regression test)", {
  set.seed(42)
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  set.seed(42)
  cond <- conditional(spmod, newdata = newexdata[1:2, ], samples = 4)

  expect_equal(
    round(as.vector(cond), 4),
    c(2.164, -0.8827, 1.665, 0.8597, -0.4249, -1.1588, 0.5448, -0.3535)
  )
})

test_that("conditional() simulate_covparams = TRUE works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_false(is.null(spmod$vcov$cov)) # ddf = "satterthwaite" default for n <= 500

  cond1 <- conditional(spmod, newdata = newexdata, samples = 50, simulate_covparams = TRUE)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))

  cond_all <- conditional(spmod, newdata = newexdata, output = "all", samples = 30, simulate_covparams = TRUE)
  expect_named(cond_all, c("newdata", "beta", "object"))
  expect_equal(dim(cond_all$newdata), c(NROW(newexdata), 30))
  expect_equal(dim(cond_all$beta), c(length(coef(spmod)), 30))
  expect_equal(dim(cond_all$object), c(spmod$n, 30))

  # samples defaults to 500 (not 10,000) under simulate_covparams = TRUE
  expect_equal(ncol(conditional(spmod, newdata = newexdata, simulate_covparams = TRUE)), 1000)
  expect_equal(ncol(conditional(spmod, newdata = newexdata)), 1000)

  # propagating covariance parameter uncertainty should not shrink the
  # marginal variance of the conditional draws relative to holding covariance
  # parameters fixed
  R <- 6000
  set.seed(1)
  cond_cp <- conditional(spmod, newdata = newexdata, samples = R, simulate_covparams = TRUE)
  set.seed(2)
  cond_fixed <- conditional(spmod, newdata = newexdata, samples = R, simulate_covparams = FALSE)
  expect_true(all(apply(cond_cp, 1, sd) > apply(cond_fixed, 1, sd)))

  # random effects and anisotropy both work through the per-draw covmatrix()
  # reuse (object_b's coefficients substituted in)
  exdata_re <- exdata
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  newexdata_re <- newexdata
  newexdata_re$grp <- factor(sample(letters[1:4], NROW(newexdata_re), replace = TRUE))
  spmod_re <- splm(y ~ x, exdata_re, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~grp)
  cond_re <- conditional(spmod_re, newdata = newexdata_re, samples = 20, simulate_covparams = TRUE)
  expect_true(all(is.finite(cond_re)))

  spmod_anis <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, anisotropy = TRUE)
  cond_anis <- conditional(spmod_anis, newdata = newexdata, samples = 20, simulate_covparams = TRUE)
  expect_true(all(is.finite(cond_anis)))

  # the exact path leaves simulate_covparams alone (no message)
  expect_no_message(
    conditional(spmod, newdata = newexdata, samples = 20, simulate_covparams = TRUE)
  )

  # a strongly-worded warning is issued above n = 500 (faking n avoids fitting
  # an actually-large model just for this check)
  spmod_fake <- spmod
  spmod_fake$n <- 600
  expect_warning(
    conditional(spmod_fake, newdata = newexdata, samples = 2, simulate_covparams = TRUE),
    "exceedingly long"
  )
  expect_no_warning(
    conditional(spmod, newdata = newexdata, samples = 2, simulate_covparams = TRUE)
  )

  # requires object$vcov$cov (i.e. ddf = "satterthwaite") -- a clear error,
  # not a cryptic one, when that covariance matrix was never computed
  spmod_asymp <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, ddf = "asymptotic")
  expect_error(
    conditional(spmod_asymp, newdata = newexdata, samples = 2, simulate_covparams = TRUE),
    "requires object\\$vcov\\$cov"
  )

  # simulate_covparams must be a single logical
  expect_error(conditional(spmod, newdata = newexdata, simulate_covparams = "yes"), "simulate_covparams must be")
})

test_that("conditional() output = 'cov'/'spcov'/'randcov' works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  spcov_names <- names(coef(spmod, type = "spcov"))

  cov_out <- conditional(spmod, newdata = newexdata, samples = 30, simulate_covparams = TRUE, output = "cov")
  expect_equal(dim(cov_out), c(length(spcov_names), 30))
  expect_equal(rownames(cov_out), spcov_names)
  expect_true(all(is.finite(cov_out)))

  # no random effects in this model -- "spcov" and "cov" coincide, "randcov" is NULL
  spcov_out <- conditional(spmod, newdata = newexdata, samples = 30, simulate_covparams = TRUE, output = "spcov")
  expect_equal(dim(spcov_out), c(length(spcov_names), 30))
  randcov_out <- conditional(spmod, newdata = newexdata, samples = 30, simulate_covparams = TRUE, output = "randcov")
  expect_null(randcov_out)

  # combining with the pre-existing outputs still works, in any combination
  combo <- conditional(spmod, newdata = newexdata,
    samples = 30, simulate_covparams = TRUE, output = c("newdata", "beta", "cov", "spcov")
  )
  expect_named(combo, c("newdata", "beta", "cov", "spcov"))
  expect_equal(dim(combo$newdata), c(NROW(newexdata), 30))
  expect_equal(dim(combo$beta), c(length(coef(spmod)), 30))
  expect_equal(dim(combo$cov), c(length(spcov_names), 30))
  expect_equal(dim(combo$spcov), c(length(spcov_names), 30))

  # requesting these without simulate_covparams = TRUE errors clearly
  expect_error(
    conditional(spmod, newdata = newexdata, samples = 5, output = "cov"),
    "only include \"cov\", \"spcov\", or \"randcov\" when simulate_covparams = TRUE"
  )

  # random effects: "randcov" is populated, and "cov" stacks spcov then randcov rows
  exdata_re <- exdata
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  newexdata_re <- newexdata
  newexdata_re$grp <- factor(sample(letters[1:4], NROW(newexdata_re), replace = TRUE))
  spmod_re <- splm(y ~ x, exdata_re, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~grp)
  randcov_names <- names(coef(spmod_re, type = "randcov"))

  randcov_out2 <- conditional(spmod_re, newdata = newexdata_re, samples = 30, simulate_covparams = TRUE, output = "randcov")
  expect_equal(dim(randcov_out2), c(length(randcov_names), 30))
  expect_equal(rownames(randcov_out2), randcov_names)

  cov_out2 <- conditional(spmod_re, newdata = newexdata_re, samples = 30, simulate_covparams = TRUE, output = "cov")
  expect_equal(rownames(cov_out2), c(names(coef(spmod_re, type = "spcov")), randcov_names))
})

test_that("simulate_theta_draw()'s clamp-to-boundary fallback is always valid", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  spcov_params_full <- coef(spmod, type = "spcov")
  spcov_type <- class(spcov_params_full)
  spcov_names_free <- names(spmod$is_known$spcov)[!spmod$is_known$spcov]
  theta_hat_free <- spcov_params_full[spcov_names_free]

  # a huge covariance matrix forces essentially every reject-and-redraw
  # attempt to fail, exercising the clamp-to-boundary fallback
  huge_vcov <- diag(1e6, length(theta_hat_free))
  dimnames(huge_vcov) <- list(spcov_names_free, spcov_names_free)
  huge_lowchol <- t(chol(huge_vcov))

  set.seed(42)
  draws <- lapply(1:30, function(i) {
    simulate_theta_draw(
      theta_hat_free, huge_lowchol, spcov_type, spcov_names_free,
      character(0), spcov_params_full, NULL,
      max_attempts = 5
    )
  })
  de_vals <- vapply(draws, function(x) x$spcov_params_b[["de"]], numeric(1))
  ie_vals <- vapply(draws, function(x) x$spcov_params_b[["ie"]], numeric(1))
  range_vals <- vapply(draws, function(x) x$spcov_params_b[["range"]], numeric(1))

  expect_true(all(is.finite(c(de_vals, ie_vals, range_vals))))
  expect_true(all(de_vals >= 0) && all(ie_vals >= 0) && all(range_vals >= 0))
  # a clamped range of exactly 0 must also force de to 0 (avoids
  # division-by-zero in exp(-d/range)-style formulas downstream)
  expect_true(all(de_vals[range_vals == 0] == 0))
})

