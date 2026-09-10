skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

# SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

# add variables
n <- NROW(exdata)
exdata$bern <- rbinom(n, size = 1, prob = 0.5)
exdata$bernfac <- factor(ifelse(exdata$bern == 0, "a", "b"))
exdata$size <- 5
exdata$bin <- rbinom(n, size = exdata$size, prob = 0.5)
exdata$prop <- runif(n)
exdata$count <- rpois(n, lambda = 5)
exdata$cont <- rgamma(n, shape = 1, rate = 1)
exdata$offset <- 1.2

test_that("the model runs for binomial data", {
  spgmod <- spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  expect_s3_class(spgmod, "spglm")
  expect_vector(AUROC(spgmod))
  spgmod <- spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  expect_s3_class(spgmod, "spglm")
  expect_vector(AUROC(spgmod))
  expect_error(spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
  spgmod <- suppressWarnings(spglm(cbind(bin, size) ~ x, family = "binomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "ml")) # ml boundary warning
  expect_s3_class(spgmod, "spglm")
  expect_error(AUROC(spgmod))
  spgmod <- spglm(y > 0 ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  expect_s3_class(spgmod, "spglm")
  expect_vector(AUROC(spgmod))


  # complicated models
  # optim non-convergence warning expected for these numerically marginal fixtures
  expect_error(suppressWarnings(spglm(bern ~ x + offset(offset), family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE)), NA)
  # check when data have an unequal number of local observations
  expect_error(spglm(bern ~ x + offset(offset), family = binomial, data = exdata[-1, , drop = FALSE], xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE), NA)
  expect_error(suppressWarnings(spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE)), NA)
  expect_error(suppressWarnings(spglm(cbind(bin, size) ~ x, family = "binomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "ml", partition_factor = ~group)), NA) # ml boundary warning

  # check glm
  spgmod1 <- spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none")
  spgmod2 <- spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential")
  expect_false(isTRUE(all.equal(as.vector(coef(spgmod1)), as.vector(coef(spgmod2)), tolerance = 0.0001)))
  gmod1 <- stats::glm(bern ~ x, family = binomial, data = exdata)
  expect_true(isTRUE(all.equal(as.vector(coef(spgmod1)), as.vector(coef(gmod1)), tolerance = 0.0001)))
})

test_that("the model runs for proportion data", {
  expect_error(spglm(prop ~ x, family = "beta", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
  expect_error(spglm(prop ~ x, family = beta, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"), NA)

  # complicated models
  expect_error(spglm(prop ~ x, family = "beta", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", anisotropy = TRUE), NA)
  expect_error(spglm(prop ~ x + offset(offset),
    family = beta, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "ie", estmethod = "reml",
    local = list(method = "kmeans")
  ), NA)
})

test_that("the model runs for count data", {
  spgmod <- spglm(count ~ x, family = poisson, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml")
  expect_s3_class(spgmod, "spglm")
  expect_error(AUROC(spgmod))
  expect_error(suppressWarnings(spglm(count ~ x, family = "nbinomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml")), NA) # ml boundary warning

  # complicated models
  expect_error(spglm(count ~ x + offset(offset),
    family = poisson, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml",
    partition_factor = ~group, local = TRUE
  ), NA)
  expect_error(suppressWarnings(spglm(count ~ x, family = "nbinomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml", random = ~group)), NA) # ml boundary warning
})

test_that("the model runs for continuous data", {
  expect_error(spglm(cont ~ x, family = "Gamma", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml"), NA)
  expect_error(spglm(cont ~ x, family = inverse.gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"), NA)

  # SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
  expect_error(spglm(cont ~ x, family = gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"))


  # complicated models
  expect_error(spglm(cont ~ x + offset(offset),
    family = "Gamma", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml",
    random = ~ group + subgroup, local = TRUE
  ), NA)
  expect_error(spglm(cont ~ x, family = inverse.gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml", partition_factor = ~subgroup, anisotropy = TRUE), NA)
  ## SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
  expect_error(spglm(cont ~ x,
    family = gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml",
    random = ~group
  ))
})

test_that("the model runs on other data sets", {
  expect_error(spglm(abs(y) ~ x, family = "Gamma", data = exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
  expect_error(suppressWarnings(spglm(abs(y) ~ x, family = Gamma, data = exdata_poly, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml")), NA) # POINT, xcoord, ycoord warnings


  # complicated models
  expect_error(spglm(abs(y) ~ x, family = "Gamma", data = exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, partition_factor = ~group, anisotropy = TRUE, local = TRUE), NA)
  expect_error(suppressWarnings(spglm(abs(y) ~ x,
    family = Gamma, data = exdata_poly, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml",
    random = ~group, partition_factor = ~subgroup
  )), NA) # POINT, xcoord, ycoord warnings
})


test_that("emmeans works", {
  spcov_type <- "exponential"
  spgmod <- spglm(abs(y) ~ x * group, family = "Gamma", exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(as.matrix(model.frame(delete.response(terms(spgmod)), spgmod$obdata)), as.matrix(emmeans::recover_data(spgmod)))
  expect_error(emmeans::emmeans(spgmod, ~ group, by = "x"), NA)
})

test_that("emmeans works missing", {
  spcov_type <- "exponential"
  spgmod <- spglm(abs(y) ~ x * group, family = "Gamma", exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(as.matrix(model.frame(delete.response(terms(spgmod)), spgmod$obdata)), as.matrix(emmeans::recover_data(spgmod)))
  expect_error(emmeans::emmeans(spgmod, ~ group, by = "x"), NA)
})

test_that("range_constrain works", {
  spcov_type <- "exponential"
  expect_error(spglm(abs(y) ~ x * group, family = "Gamma", exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", range_constrain = TRUE), NA)
  expect_error(suppressWarnings(spglm(abs(y) ~ x * group, family = "Gamma", exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", range_constrain = TRUE)), NA) # ml boundary warning
})

test_that("ml boundary regression: an optimizer that wanders to the floor is reconciled and flagged", {
  # reproduces the scenario documented in the ml-boundary write-up: data simulated
  # with genuine spatial structure, whose exponential ml fit still converges to a
  # tiny ie chasing the boundary artifact
  set.seed(42)
  n_side <- 12
  coords <- expand.grid(xcoord = seq_len(n_side), ycoord = seq_len(n_side))
  coords$x <- rnorm(nrow(coords))

  sp_params <- spcov_params("exponential", de = 0.6, ie = 0.1, range = 4, rotate = 0, scale = 1)
  coords$y <- sprpois(sp_params, mean = 0.5 + 0.3 * coords$x, data = coords, xcoord = xcoord, ycoord = ycoord)

  expect_warning(
    mod_exp <- spglm(y ~ x, family = "poisson", data = coords, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"),
    NA # de stays substantial (see the write-up), so this specific fit should not trip the boundary warning
  )
  # the optimizer wanders below the floor for ie; the reported value must be reconciled to it
  expect_true(coef(mod_exp, type = "spcov")[["ie"]] >= 1e-4)

  expect_warning(
    mod_none <- spglm(y ~ x, family = "poisson", data = coords, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "ml"),
    "numerical boundary"
  )
  # reml should not show the same AIC preference for "none" over "exponential"
  mod_exp_reml <- spglm(y ~ x, family = "poisson", data = coords, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  mod_none_reml <- spglm(y ~ x, family = "poisson", data = coords, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml")
  expect_true(AIC(mod_exp_reml) < AIC(mod_none_reml))
})

test_that("fitted-probability saturation warning fires under spatial separation but not on a sane binomial fit", {
  skip_if_not_installed("sf")

  # reproduces the scenario documented in the separation write-up: default
  # (reml) estmethod, de runs away, fitted probabilities fully saturate
  data(lake)
  lake_df <- sf::st_drop_geometry(lake)
  xy <- sf::st_coordinates(lake)
  lake_df$xcoord <- xy[, 1]
  lake_df$ycoord <- xy[, 2]
  lake_df$y <- as.numeric(lake_df$log_cond > 3)

  # this fixture also fails to converge (the saturation and non-convergence
  # diagnostics commonly co-occur); nest so both expected warnings are caught
  expect_warning(
    expect_warning(
      spglm(y ~ temp + elev, data = lake_df, family = binomial, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential"),
      "did not converge"
    ),
    "separation"
  )

  # a sane binomial fit (well below the saturation threshold, see the
  # calibration in the separation write-up) must not trip the warning
  expect_warning(
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"),
    NA
  )
})

test_that("ml boundary warning fires for spcov_type = none but not reml", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "ml"),
    "numerical boundary"
  )
  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"),
    NA
  )
})

test_that("ml boundary warning uses the actual fitted value and does not override a known ie", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  spinit <- spcov_initial("exponential", de = 1e-8, ie = 1e-8, range = 1, known = c("de", "ie", "range"))
  expect_warning(
    spmod <- spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spinit, estmethod = "ml"),
    "numerical boundary"
  )
  # a known (fixed) ie must never be reconciled to the diagtol floor
  expect_equal(coef(spmod, type = "spcov")[["ie"]], 1e-8)
})

test_that("no ml boundary warning for a well-identified spatial covariance", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"),
    NA
  )
})

test_that("optim non-convergence warning fires for spglm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", control = list(maxit = 1)),
    "did not converge"
  )
  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential"),
    NA
  )
})

test_that("anisotropy = TRUE with fully known covariance and dispersion parameters does not error", {
  # regression test for a bug where run_laploglik_dispatch_spglm() called
  # use_laploglik_known_anis() with an unnamed positional argument list one
  # short of its formals, silently dropping randcov_initial and erroring with
  # "argument \"randcov_initial\" is missing, with no default" -- for both the
  # no-random-effects case and the with-random-effects case
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  spcov_initial_known <- spcov_initial("exponential",
                                       de = 1, ie = 0.5, range = 2, rotate = 0.3, scale = 0.7,
                                       known = "given"
  )
  dispersion_initial_known <- dispersion_initial("poisson", dispersion = 1, known = "given")

  # no random effects
  expect_error(
    spmod_no_randcov <- spglm(y_pois ~ x,
                              family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord,
                              spcov_initial = spcov_initial_known,
                              dispersion_initial = dispersion_initial_known,
                              anisotropy = TRUE
    ),
    NA
  )
  expect_equal(as.numeric(coef(spmod_no_randcov, type = "spcov")), as.numeric(spcov_initial_known$initial))

  # with a known random effect
  randcov_initial_known <- randcov_initial(group = 0.2, known = "given")
  expect_error(
    spmod_with_randcov <- spglm(y_pois ~ x,
                                family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord,
                                spcov_initial = spcov_initial_known,
                                dispersion_initial = dispersion_initial_known,
                                random = ~group,
                                randcov_initial = randcov_initial_known,
                                anisotropy = TRUE
    ),
    NA
  )
  expect_equal(as.numeric(coef(spmod_with_randcov, type = "spcov")), as.numeric(spcov_initial_known$initial))
  expect_equal(unname(coef(spmod_with_randcov, type = "randcov")[["1 | group"]]), 0.2)
})

test_that("predict() errors informatively when newdata is missing coordinate columns", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod1 <- spglm(abs(y) ~ x, "Gamma", exdata, spcov_type = "exponential", xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml")

  # sanity: predict() still works normally when both coordinate columns are present
  expect_vector(predict(spmod1, newdata = newexdata))

  newexdata_nocoord <- newexdata
  newexdata_nocoord$xcoord <- NULL
  newexdata_nocoord$ycoord <- NULL
  expect_error(predict(spmod1, newdata = newexdata_nocoord), "coordinate column")
})

test_that("predict() errors informatively when newdata has NA in a random intercept grouping column", {
  set.seed(11)
  n <- 60
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  grp <- factor(sample(letters[1:5], n, replace = TRUE))
  p <- plogis(0.3 * x1)
  ybin <- rbinom(n, 1, p)
  d <- data.frame(ybin = ybin, x1 = x1, grp = grp, xcoord = xcoord, ycoord = ycoord)

  gmod <- spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~ (1 | grp))

  newdata_good <- data.frame(
    x1 = c(0.5, -0.5), grp = factor(c("a", "b"), levels = levels(grp)),
    xcoord = c(0.3, 0.6), ycoord = c(0.3, 0.6)
  )
  expect_vector(predict(gmod, newdata_good))

  newdata_na <- newdata_good
  newdata_na$grp[1] <- NA
  expect_error(predict(gmod, newdata_na), "Cannot have NA values in predictors.")
  expect_error(predict(gmod, newdata_na, se.fit = TRUE), "Cannot have NA values in predictors.")
})

test_that("spglm() errors informatively when a formula/random/partition_factor variable is not in data", {
  set.seed(12)
  n <- 30
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  grp <- factor(sample(letters[1:3], n, replace = TRUE))
  ybin <- rbinom(n, 1, plogis(0.3 * x1))
  d <- data.frame(ybin = ybin, x1 = x1, grp = grp, xcoord = xcoord, ycoord = ycoord)

  # a same-named object in the calling environment (but not in data) should
  # not be silently picked up via ordinary formula scoping -- it should error
  not_a_col <- rnorm(n)
  not_a_group <- factor(sample(letters[1:3], n, replace = TRUE))

  expect_error(spglm(ybin ~ not_a_col, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "not_a_col.*not found in data")
  expect_error(spglm(not_a_col ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "not_a_col.*not found in data")
  expect_error(spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~not_a_group), "not_a_group.*not found in data")
  expect_error(spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, partition_factor = ~not_a_group), "not_a_group.*not found in data")

  # sanity: a valid call still works
  expect_s3_class(spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "spglm")
})

test_that("spglm() formula supports . as shorthand for all predictors, excluding coordinates", {
  # see the matching test in test-splm.R for the full rationale
  set.seed(13)
  n <- 30
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  ybin <- rbinom(n, 1, plogis(0.3 * x1))
  d <- data.frame(ybin = ybin, x1 = x1, xcoord = xcoord, ycoord = ycoord)

  mod_dot <- spglm(ybin ~ ., family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  mod_explicit <- spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_equal(names(coef(mod_dot)), names(coef(mod_explicit)))
  expect_equal(unname(coef(mod_dot)), unname(coef(mod_explicit)))
  expect_false(any(c("xcoord", "ycoord") %in% colnames(model.matrix(mod_dot))))

  d_sf <- sf::st_as_sf(d, coords = c("xcoord", "ycoord"))
  mod_sf_dot <- spglm(ybin ~ ., family = "binomial", data = d_sf, spcov_type = "exponential")
  expect_equal(unname(coef(mod_sf_dot)), unname(coef(mod_explicit)))

  expect_error(
    spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~.),
    "not supported in random"
  )
})

test_that("a user-supplied local$index must match the length of the non-missing response", {
  # see the matching test in test-splm.R for the full rationale
  set.seed(12)
  n <- 40
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  ybin <- rbinom(n, 1, plogis(0.3 * x1))
  d <- data.frame(ybin = ybin, x1 = x1, xcoord = xcoord, ycoord = ycoord)

  na_rows <- c(3, 10, 17, 25, 33)
  d_na <- d
  d_na$ybin[na_rows] <- NA
  n_obs <- n - length(na_rows)

  full_index <- sample(1:4, n, replace = TRUE) # wrong length: matches original n, not n_obs
  observed_index <- which(!is.na(d_na$ybin))
  correct_index <- full_index[observed_index] # right length: matches n_obs

  expect_error(
    get_data_object_spglm(
      formula = ybin ~ x1, family = "binomial", data = d_na, spcov_initial = spcov_initial("exponential"),
      xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml", anisotropy = FALSE,
      random = NULL, randcov_initial = NULL, partition_factor = NULL,
      local = list(index = full_index), range_constrain = FALSE
    ),
    "local\\$index must have the same length"
  )
  expect_error(
    spglm(ybin ~ x1, family = "binomial", data = d_na, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = list(index = full_index)),
    "local\\$index must have the same length"
  )

  data_object <- get_data_object_spglm(
    formula = ybin ~ x1, family = "binomial", data = d_na, spcov_initial = spcov_initial("exponential"),
    xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml", anisotropy = FALSE,
    random = NULL, randcov_initial = NULL, partition_factor = NULL,
    local = list(index = correct_index), range_constrain = FALSE
  )
  expect_length(data_object$local_index, n_obs)
  expect_equal(as.vector(data_object$local_index), as.vector(correct_index))

  expect_s3_class(
    spglm(ybin ~ x1, family = "binomial", data = d_na, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = list(index = correct_index)),
    "spglm"
  )
})
