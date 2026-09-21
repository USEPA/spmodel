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

test_that("the model runs for exponential", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for exponential (partition group)", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", partition_factor = ~group), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", partition_factor = ~group), NA)
})

test_that("the model runs for exponential (random group)", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", random = ~group), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~group), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  randcov_initial_val <- randcov_initial(group = 1)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~group, randcov_initial = randcov_initial_val), NA)
})

test_that("the model runs for exponential (random and subgroup)", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~ group + subgroup), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", random = ~ group + subgroup), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~ group + subgroup), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~group), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  randcov_initial_val <- randcov_initial(group = 1)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~ group + subgroup, randcov_initial = randcov_initial_val), NA)
})

test_that("the model runs for exponential (random nested subgroup)", {
  spcov_type <- "exponential"
  # a nested subgroup random effect leaves little information to estimate
  # its own covariance parameter, which regularly makes the automatic
  # (n <= 500) Satterthwaite ddf computation's covariance-parameter Hessian
  # non-positive-definite -- spmodel already warns and falls back to NULL
  # ddf gracefully, so the warning is expected/tolerated. The Hessian's
  # positive-definiteness is borderline for every call in this block (not
  # just some), so whether a given call warns can flip between runs/
  # machines (BLAS-level floating-point nondeterminism) -- suppress across
  # the board rather than pin to whichever calls happened to warn once
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~ group / subgroup)), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", random = ~ group / subgroup)), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~ group / subgroup)), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~ group / subgroup)), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  randcov_initial_val <- randcov_initial(group = 1)
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~ group / subgroup, randcov_initial = randcov_initial_val)), NA)
})

test_that("nested random effect (group/subgroup) uses the full crossed grouping, not just group (regression test)", {
  # bug test where the "group:subgroup" term's design matrix
  # should not collapse to only the first group's levels (instead of group/subgroup combinations)
  n_combos_observed <- nlevels(droplevels(interaction(exdata$group, exdata$subgroup)))
  expect_true(n_combos_observed > nlevels(exdata$group))

  Z_nested <- get_randcov_Z("1 | group:subgroup", exdata)$Z
  expect_equal(ncol(Z_nested), n_combos_observed)

  spmod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", random = ~ group / subgroup)
  expect_true("1 | group:subgroup" %in% names(coef(spmod, type = "randcov")))

  # fitted(type = "randcov") returns one BLUP per level of each random
  # effect term
  nested_blups <- fitted(spmod, type = "randcov")[["1 | group:subgroup"]]
  expect_equal(length(nested_blups), n_combos_observed)
  expect_equal(sort(names(nested_blups)), sort(colnames(Z_nested)))
})

test_that("the model runs for exponential (random and partitioning)", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group, partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", random = ~group, partition_factor = ~group), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group, partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~group, partition_factor = ~group), NA)
})

test_that("the model runs for anisotropy", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", anisotropy = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", anisotropy = TRUE
  ), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  # de known at a fixed value leaves the automatic Satterthwaite ddf
  # computation with less information than usual about the remaining
  # covariance parameters, which regularly makes its covariance-parameter
  # Hessian non-positive-definite -- spmodel already warns and falls back to
  # NULL ddf gracefully, so the warning is expected/tolerated here
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "reml", anisotropy = TRUE
  )), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "ml", anisotropy = TRUE
  )), NA)
})

test_that("the model runs for and random effects", {
  spcov_type <- "exponential"
  # anisotropy + a random effect leaves the automatic Satterthwaite ddf
  # computation with a lot of covariance parameters to estimate relative to
  # this small fixture, which regularly makes its covariance-parameter
  # Hessian non-positive-definite -- spmodel already warns and falls back to
  # NULL ddf gracefully, so the warning is expected/tolerated for every call
  # in this block
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", anisotropy = TRUE, random = ~group
  )), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", anisotropy = TRUE, random = ~group
  )), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "reml", anisotropy = TRUE, random = ~group
  )), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "ml", anisotropy = TRUE, random = ~group
  )), NA)
})

test_that("the model runs for and partitioning", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", anisotropy = TRUE, partition_factor = ~group
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", anisotropy = TRUE, partition_factor = ~group
  ), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  # de known at a fixed value leaves the automatic Satterthwaite ddf
  # computation with less information than usual about the remaining
  # covariance parameters, which regularly makes its covariance-parameter
  # Hessian non-positive-definite -- spmodel already warns and falls back to
  # NULL ddf gracefully, so the warning is expected/tolerated here
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "reml", anisotropy = TRUE, partition_factor = ~group
  )), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "ml", anisotropy = TRUE, partition_factor = ~group
  )), NA)
})

test_that("the model runs for and random effects and partitioning", {
  spcov_type <- "exponential"
  # optim non-convergence warning expected for this numerically marginal fixture
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", anisotropy = TRUE, random = ~group, partition_factor = ~group
  )), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", anisotropy = TRUE, random = ~group, partition_factor = ~group
  ), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  # de known plus anisotropy, a random effect, and partitioning leaves the
  # automatic Satterthwaite ddf computation with a lot of covariance
  # parameters to estimate relative to this small fixture, which regularly
  # makes its covariance-parameter Hessian non-positive-definite -- spmodel
  # already warns and falls back to NULL ddf gracefully, so the warning is
  # expected/tolerated here
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "reml",
    anisotropy = TRUE, random = ~group, partition_factor = ~group
  )), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "ml",
    anisotropy = TRUE, random = ~group, partition_factor = ~group
  )), NA)
})

test_that("the model runs for exponential and missing data", {
  spcov_type <- spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})


test_that("the model runs for exponential (random group) and missing data", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group), NA)
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", random = ~group), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group), NA)
  expect_error(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~group), NA)
})


test_that("the model runs for big data", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", local = TRUE), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", local = TRUE), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", local = TRUE), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl", local = TRUE), NA)

  # parallel for REML and ML and no errors for other methods
  # CRAN ONLY ALLOWS 2 CORES FOR TESTING
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", local = list(parallel = TRUE, ncores = 2)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", local = list(parallel = TRUE, ncores = 2)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-wls", local = list(parallel = TRUE, ncores = 2)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-cl", local = list(parallel = TRUE, ncores = 2)
  ), NA)

  # in case var_adjust default changed to "none"
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", local = list(parallel = TRUE, var_adjust = "none", ncores = 2)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", local = list(parallel = TRUE, var_adjust = "theoretical", ncores = 2)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-wls", local = list(parallel = TRUE, var_adjust = "pooled", ncores = 2)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-cl", local = list(parallel = TRUE, var_adjust = "empirical", ncores = 2)
  ), NA)

  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", local = list(size = 30)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", local = list(size = 30)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-wls", local = list(size = 30)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-cl", local = list(size = 30)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", local = list(groups = 10)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", local = list(groups = 10)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-wls", local = list(groups = 10)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-cl", local = list(groups = 10)
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", local = list(index = sample(1:4, size = 100, replace = TRUE))
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", local = list(index = sample(1:4, size = 100, replace = TRUE))
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-wls", local = list(index = sample(1:4, size = 100, replace = TRUE))
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-cl", local = list(index = sample(1:4, size = 100, replace = TRUE))
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", local = list(method = "kmeans")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", local = list(method = "kmeans")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-wls", local = list(method = "kmeans")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-cl", local = list(method = "kmeans")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", local = list(var_adjust = "none")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", local = list(var_adjust = "theoretical")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-wls", local = list(var_adjust = "empirical")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-cl", local = list(var_adjust = "pooled")
  ), NA)


  # random effects
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", random = ~group, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", random = ~group, local = TRUE
  ), NA)

  # random effects (nested)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", random = ~ group / subgroup, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", random = ~ group / subgroup, local = TRUE
  ), NA)


  # random effects (x2)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", random = ~ group + subgroup, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", random = ~ group + subgroup, local = TRUE
  ), NA)

  # anisotropy
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", anisotropy = TRUE, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", anisotropy = TRUE, local = TRUE
  ), NA)


  # partitioning
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", partition_factor = ~group, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", partition_factor = ~group, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-wls", partition_factor = ~group, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "sv-cl", partition_factor = ~group, local = TRUE
  ), NA)

  # random effects partitioning
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", random = ~group, partition_factor = ~group, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", random = ~group, partition_factor = ~group, local = TRUE
  ), NA)



  # random effects anisotropy
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", random = ~group, anisotropy = TRUE, local = TRUE
  ), NA)

  # partitioning anisotropy
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", partition_factor = ~group, anisotropy = TRUE, local = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", partition_factor = ~group, anisotropy = TRUE, local = TRUE
  ), NA)

  # random effects partitioning anisotropy
  # optim non-convergence warning expected for this numerically marginal fixture
  expect_error(suppressWarnings(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml",
    random = ~group, partition_factor = ~group, anisotropy = TRUE, local = TRUE
  )), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml",
    random = ~group, partition_factor = ~group, anisotropy = TRUE, local = TRUE
  ), NA)
})

test_that("the model runs for spherical", {
  spcov_type <- "spherical"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for gaussian", {
  spcov_type <- "gaussian"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for triangular", {
  spcov_type <- "triangular"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)

  # try giving y coordinate
  # ddf = "asymptotic": these calls are already under test for a specific,
  # unrelated warning (providing ycoord for a 1D-only covariance);
  # expect_warning() only consumes the first warning it sees and lets any
  # other warning from the same call leak through unmuffled, and this small
  # fixture's automatic (n <= 500) Satterthwaite ddf computation regularly
  # emits a second, incidental "not positive definite" warning that would
  # otherwise leak through that way -- suppressWarnings() isn't an option
  # here since it would also swallow the warning expect_warning() needs to
  # see, so ddf is disabled directly instead
  expect_warning(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", ddf = "asymptotic"))
  expect_warning(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", ddf = "asymptotic"))
})

test_that("the model runs for circular", {
  spcov_type <- "circular"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for none", {
  spcov_type <- "none"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for ie", {
  spcov_type <- "ie"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
  mod1 <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
  mod2 <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml")
  expect_equal(as.vector(coef(mod1, "spcov")), as.vector(coef(mod2, "spcov")), tolerance = 0.01)
})

test_that("known Gaussian iid variance is retained for none and ie", {
  dat <- exdata
  dat$off <- 0.2 * sin(seq_len(nrow(dat)))
  X <- model.matrix(y ~ x + offset(off), dat)
  expected_vcov <- 0.2 * solve(crossprod(X))

  for (estmethod in c("ml", "reml")) {
    for (spcov_type in c("none", "ie")) {
      fit <- splm(y ~ x + offset(off), dat,
        estmethod = estmethod, ddf = "asymptotic",
        spcov_initial = spcov_initial(spcov_type, ie = 0.2, known = "given")
      )
      expect_equal(coef(fit, "spcov")[["ie"]], 0.2)
      expect_equal(as.matrix(covmatrix(fit)), diag(0.2, nrow(dat)))
      expect_equal(vcov(fit), expected_vcov, tolerance = 1e-10)
    }
  }
})

test_that("numerical nugget floor is excluded from BLUPs and included in prediction variance", {
  n_obs <- 45
  dat <- data.frame(
    xc = seq(0, by = 1000, length.out = n_obs),
    yc = 700 * sin(seq_len(n_obs) / 4),
    x = seq(-1, 1, length.out = n_obs),
    off = 0.2 * cos(seq_len(n_obs) / 5)
  )
  dat$y <- 1 + 0.7 * dat$x + dat$off + 0.35 * sin(seq_len(n_obs) / 3)
  newdata <- data.frame(
    xc = c(2500, 8500, 17500, 32500),
    yc = c(300, -450, 600, -200),
    x = c(-0.8, -0.25, 0.35, 0.8),
    off = c(0.1, -0.05, 0.15, -0.1)
  )
  de <- 2
  range <- 12000
  dist_obs <- sqrt(
    outer(dat$xc, dat$xc, "-")^2 + outer(dat$yc, dat$yc, "-")^2
  )
  dist_pred <- sqrt(
    outer(newdata$xc, dat$xc, "-")^2 + outer(newdata$yc, dat$yc, "-")^2
  )
  K <- de * exp(-dist_obs / range)
  C <- de * exp(-dist_pred / range)

  for (ie in c(0, 1e-9, 0.2)) {
    fit <- splm(y ~ x + offset(off), dat,
      xcoord = xc, ycoord = yc, estmethod = "reml", ddf = "asymptotic",
      spcov_initial = spcov_initial("exponential",
        de = de, ie = ie, range = range, known = "given"
      )
    )
    effective_ie <- max(ie, 1e-4 * de)
    V <- K + diag(effective_ie, n_obs)
    V_inv <- solve(V)
    X <- model.matrix(fit)
    y <- model.response(model.frame(fit)) - model.offset(model.frame(fit))
    B <- solve(crossprod(X, V_inv %*% X))
    beta_reference <- B %*% crossprod(X, V_inv %*% y)
    residual_weight <- V_inv %*% (y - X %*% beta_reference)

    expect_equal(coef(fit, "spcov")[["ie"]], ie)
    expect_equal(unname(coef(fit)), as.numeric(beta_reference), tolerance = 1e-10)
    expect_equal(
      unname(fitted(fit, "spcov")$de),
      as.numeric(K %*% residual_weight),
      tolerance = 1e-10
    )
    expect_equal(
      unname(fitted(fit, "spcov")$ie),
      as.numeric(ie * residual_weight),
      tolerance = 1e-10
    )

    X0 <- model.matrix(delete.response(terms(fit)), newdata)
    H <- X0 - C %*% V_inv %*% X
    fit_reference <- as.numeric(
      X0 %*% beta_reference + C %*% residual_weight + newdata$off
    )
    var_reference <- de + effective_ie -
      rowSums((C %*% V_inv) * C) + rowSums((H %*% B) * H)
    full_prediction <- predict(fit, newdata, se.fit = TRUE, local = FALSE)
    expect_equal(unname(full_prediction$fit), fit_reference, tolerance = 1e-10)
    expect_equal(unname(full_prediction$se.fit), sqrt(as.numeric(var_reference)), tolerance = 1e-10)

    local_prediction <- predict(fit, newdata, se.fit = TRUE,
      local = list(method = "distance", size = 12, parallel = FALSE)
    )
    local_reference <- lapply(seq_len(nrow(newdata)), function(i) {
      keep <- order(dist_pred[i, ])[seq_len(12)]
      V_local <- V[keep, keep, drop = FALSE]
      C_local <- C[i, keep, drop = FALSE]
      X_local <- X[keep, , drop = FALSE]
      y_local <- y[keep]
      V_local_inv <- solve(V_local)
      residual_local <- y_local - X_local %*% beta_reference
      H_local <- X0[i, , drop = FALSE] - C_local %*% V_local_inv %*% X_local
      list(
        fit = as.numeric(
          X0[i, , drop = FALSE] %*% beta_reference +
            C_local %*% V_local_inv %*% residual_local + newdata$off[i]
        ),
        var = as.numeric(
          de + effective_ie -
            C_local %*% V_local_inv %*% t(C_local) + H_local %*% B %*% t(H_local)
        )
      )
    })
    expect_equal(
      unname(local_prediction$fit),
      vapply(local_reference, `[[`, numeric(1), "fit"),
      tolerance = 1e-10
    )
    expect_equal(
      unname(local_prediction$se.fit),
      sqrt(vapply(local_reference, `[[`, numeric(1), "var")),
      tolerance = 1e-10
    )
  }

  params <- spcov_params("exponential", de = de, ie = 0, range = range)
  expect_equal(as.matrix(cov_matrix_cross(params, dist_pred)), C)

  dat$part <- factor(rep(c("a", "b", "c"), each = 15))
  partition_fit <- splm(y ~ x + offset(off), dat,
    xcoord = xc, ycoord = yc, partition_factor = ~part,
    estmethod = "reml", ddf = "asymptotic",
    spcov_initial = spcov_initial("exponential",
      de = de, ie = 0, range = range, known = "given"
    )
  )
  same_partition <- outer(dat$part, dat$part, "==")
  K_partition <- K * same_partition
  V_partition <- K_partition + diag(1e-4 * de, n_obs)
  X_partition <- model.matrix(partition_fit)
  y_partition <- model.response(model.frame(partition_fit)) -
    model.offset(model.frame(partition_fit))
  residual_weight_partition <- solve(
    V_partition,
    y_partition - X_partition %*% coef(partition_fit)
  )
  expect_equal(
    unname(fitted(partition_fit, "spcov")$de),
    as.numeric(K_partition %*% residual_weight_partition),
    tolerance = 1e-10
  )
})


test_that("the model runs for cubic", {
  spcov_type <- "cubic"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for pentaspherical", {
  spcov_type <- "pentaspherical"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for cosine", {
  spcov_type <- "cosine"
  # reml/ml fits (unlike sv-wls/sv-cl) attempt the automatic n <= 500
  # Satterthwaite ddf computation at fit time; this fixture regularly makes
  # its covariance-parameter Hessian non-positive-definite -- spmodel
  # already warns and falls back to NULL ddf gracefully, so the warning is
  # expected/tolerated rather than a sign these fits are broken
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, spcov_type = spcov_type, estmethod = "reml")), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, spcov_initial = spcov_initial_val, estmethod = "reml")), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, spcov_initial = spcov_initial_val, estmethod = "ml")), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)

  # try giving y coordinate
  # ddf = "asymptotic": these calls are already under test for a specific,
  # unrelated warning (providing ycoord for a 1D-only covariance);
  # expect_warning() only consumes the first warning it sees and lets any
  # other warning from the same call leak through unmuffled, and this small
  # fixture's automatic (n <= 500) Satterthwaite ddf computation regularly
  # emits a second, incidental "not positive definite" warning that would
  # otherwise leak through that way -- suppressWarnings() isn't an option
  # here since it would also swallow the warning expect_warning() needs to
  # see, so ddf is disabled directly instead
  expect_warning(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", ddf = "asymptotic"))
  expect_warning(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", ddf = "asymptotic"))
})

test_that("the model runs for wave", {
  spcov_type <- "wave"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for jbessel", {
  spcov_type <- "jbessel"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for gravity", {
  spcov_type <- "gravity"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for rquad", {
  spcov_type <- "rquad"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for magnetic", {
  spcov_type <- "magnetic"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for matern", {
  spcov_type <- "matern"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("the model runs for cauchy", {
  spcov_type <- "cauchy"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})

test_that("cauchy covariance is stable for large range and shape", {
  ordinary_params <- spcov_params("cauchy",
    de = 2.3, ie = 0.4, range = 1.7, extra = 0.8
  )
  dist_vector <- c(0, 0.1, 1, 5)
  ordinary_reference <- ordinary_params[["de"]] *
    (1 + (dist_vector / ordinary_params[["range"]])^2)^(-ordinary_params[["extra"]])
  expect_equal(spcov_vector(ordinary_params, dist_vector), ordinary_reference)

  dist_matrix <- as.matrix(dist(dist_vector))
  matrix_reference <- ordinary_params[["de"]] *
    (1 + (dist_matrix / ordinary_params[["range"]])^2)^(-ordinary_params[["extra"]])
  diag(matrix_reference) <- diag(matrix_reference) + ordinary_params[["ie"]]
  dense_result <- spcov_matrix(ordinary_params, dist_matrix)
  expect_equal(dense_result, matrix_reference)
  expect_equal(
    as.matrix(spcov_matrix(ordinary_params, Matrix::Matrix(dist_matrix, sparse = TRUE))),
    matrix_reference
  )

  extreme_params <- spcov_params("cauchy",
    de = 1, ie = 0, range = 1e9, extra = 1e16
  )
  extreme_distances <- matrix(c(0, 1, 1, 0), 2)
  gaussian_reference <- exp(-0.01)
  extreme_matrix <- spcov_matrix(extreme_params, extreme_distances)
  expect_equal(extreme_matrix[1, 2], gaussian_reference, tolerance = 1e-12)
  expect_equal(spcov_vector(extreme_params, 1), gaussian_reference, tolerance = 1e-12)
  expect_equal(spcov_vector(extreme_params, 0), 1)
  expect_equal(diag(extreme_matrix), rep(1.0001, 2))

  gaussian_range <- 10
  shape <- c(1e4, 1e8, 1e16)
  limit_distance <- 2
  limit_values <- vapply(shape, function(extra) {
    params <- spcov_params("cauchy",
      de = 1, ie = 0, range = gaussian_range * sqrt(extra), extra = extra
    )
    spcov_vector(params, limit_distance)
  }, numeric(1))
  limit_reference <- exp(-(limit_distance / gaussian_range)^2)
  expect_true(all(diff(abs(limit_values - limit_reference)) < 0))
  expect_equal(tail(limit_values, 1), limit_reference, tolerance = 1e-15)
})

test_that("known cauchy fits match the anisotropic Gaussian limit", {
  dat <- data.frame(
    y = c(1.1, 2.2, 1.4, 2.8, 1.9, 3.2, 2.4, 3.5),
    x = c(-1, -0.6, -0.2, 0.1, 0.4, 0.7, 1, 1.3),
    xc = c(0, 1, 2, 4, 5, 7, 8, 10),
    yc = c(0, 2, 1, 3, 6, 5, 9, 8)
  )
  newdata <- data.frame(
    x = c(-0.4, 0.3, 1.1),
    xc = c(1.5, 4.5, 9),
    yc = c(1, 4, 7)
  )
  cauchy_initial <- spcov_initial("cauchy",
    de = 0.8, ie = 0.2, range = 1e9, extra = 1e16,
    rotate = 0.35, scale = 0.6, known = "given"
  )
  gaussian_initial <- spcov_initial("gaussian",
    de = 0.8, ie = 0.2, range = 10,
    rotate = 0.35, scale = 0.6, known = "given"
  )
  cauchy_fit <- splm(y ~ x, dat,
    xcoord = xc, ycoord = yc, spcov_initial = cauchy_initial,
    anisotropy = TRUE, estmethod = "reml", ddf = "asymptotic"
  )
  gaussian_fit <- splm(y ~ x, dat,
    xcoord = xc, ycoord = yc, spcov_initial = gaussian_initial,
    anisotropy = TRUE, estmethod = "reml", ddf = "asymptotic"
  )

  expect_equal(coef(cauchy_fit, "spcov"), cauchy_initial$initial)
  expect_equal(logLik(cauchy_fit), logLik(gaussian_fit), tolerance = 1e-12)
  expect_equal(
    as.matrix(covmatrix(cauchy_fit)),
    as.matrix(covmatrix(gaussian_fit)),
    tolerance = 1e-12
  )
  for (cov_type in c("pred.obs", "pred.pred")) {
    expect_equal(
      as.matrix(covmatrix(cauchy_fit, newdata, cov_type = cov_type)),
      as.matrix(covmatrix(gaussian_fit, newdata, cov_type = cov_type)),
      tolerance = 1e-12
    )
  }
  expect_equal(
    predict(cauchy_fit, newdata, se.fit = TRUE),
    predict(gaussian_fit, newdata, se.fit = TRUE),
    tolerance = 1e-12
  )
})

test_that("the model runs for pexponential", {
  spcov_type <- "pexponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)
})


test_that("the model runs for all sv-wls weights", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", weights = "cressie"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", weights = "cressie-dr"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", weights = "cressie-nopairs"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", weights = "cressie-dr-nopairs"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", weights = "pairs"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", weights = "pairs-invd"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", weights = "pairs-invrd"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", weights = "ols"), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", weights = "cressie"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", weights = "cressie-dr"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", weights = "cressie-nopairs"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", weights = "cressie-dr-nopairs"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", weights = "pairs"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", weights = "pairs-invd"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", weights = "pairs-invrd"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", weights = "ols"), NA)
})


test_that("the model runs for certain known parameter assignments", {
  # 3 param geo
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl"), NA)
  randcov_initial_val <- randcov_initial(group = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml", random = ~group, randcov_initial = randcov_initial_val), NA)
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, rotate = 2, scale = 0.5, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml", random = ~group, randcov_initial = randcov_initial_val), NA)

  # 2 param geo
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 0, range = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl"), NA)
  randcov_initial_val <- randcov_initial(group = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 0, range = 1, rotate = 2, scale = 0.5, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)

  # 4 param geo
  spcov_initial_val <- spcov_initial("matern", de = 1, ie = 1, range = 1, extra = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl"), NA)
  randcov_initial_val <- randcov_initial(group = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)
  spcov_initial_val <- spcov_initial("matern", de = 1, ie = 1, range = 1, extra = 1, rotate = 2, scale = 0.5, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)

  spcov_initial_val <- spcov_initial("matern", de = 1, ie = 0, range = 1, extra = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl"), NA)
  randcov_initial_val <- randcov_initial(group = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)
  spcov_initial_val <- spcov_initial("matern", de = 1, ie = 0, range = 1, rotate = 2, extra = 1, scale = 0.5, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl"), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)

  # 3 param auto
  spcov_initial_val <- spcov_initial("car", de = 1, ie = 0, range = 0.5, extra = 1, known = "given")
  expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val), NA)
  randcov_initial_val <- randcov_initial(group = 1, known = "given")
  expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val), NA)
  spcov_initial_val <- spcov_initial("car", de = 1, ie = 0, range = 0.5, extra = 1, rotate = 2, scale = 0.5, known = "given")
  expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val), NA)
  expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val), NA)

  spcov_initial_val <- spcov_initial("car", de = 1, ie = 1, range = 0.5, extra = 1, known = "given")
  expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val), NA)
  randcov_initial_val <- randcov_initial(group = 1, known = "given")
  expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val), NA)
  spcov_initial_val <- spcov_initial("car", de = 1, ie = 1, range = 0.5, rotate = 2, extra = 1, scale = 0.5, known = "given")
  expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val), NA)
  expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val), NA)

  # 3 param geo w/ partition factors
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls", partition_factor = ~group, local = TRUE), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl", partition_factor = ~group, local = TRUE), NA)
  randcov_initial_val <- randcov_initial(group = 1, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val, partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml", random = ~group, randcov_initial = randcov_initial_val, partition_factor = ~group), NA)
  spcov_initial_val <- spcov_initial("exponential", de = 1, ie = 1, range = 1, rotate = 2, scale = 0.5, known = "given")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls", partition_factor = ~group, local = TRUE), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl", partition_factor = ~group, local = TRUE), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, random = ~group, randcov_initial = randcov_initial_val, partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord, estmethod = "ml", random = ~group, randcov_initial = randcov_initial_val, partition_factor = ~group), NA)
})


test_that("the model runs for sf and sp objects", {

  # point data
  exdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"), crs = 5070)
  exdata_sf_geo <- sf::st_transform(exdata_sf, crs = 4326)
  exdata_sf_NA <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"), crs = NA)

  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_type <- "none"
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_type <- "matern"
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_sf, spcov_type = spcov_type, estmethod = "sv-cl"), NA)

  # polygon data
  spcov_type <- "exponential"
  expect_error(suppressWarnings(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "reml")), NA) # POINT warning
  expect_error(suppressWarnings(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "ml")), NA)  # POINT warning
  expect_error(suppressWarnings(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-wls")), NA)  # POINT warning
  expect_error(suppressWarnings(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-cl")), NA)  # POINT warning
  spcov_type <- "none"
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_type <- "matern"
  expect_error(suppressWarnings(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "reml")), NA)  # POINT warning
  expect_error(suppressWarnings(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "ml")), NA)  # POINT warning
  expect_error(suppressWarnings(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-wls")), NA)  # POINT warning
  expect_error(suppressWarnings(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-cl")), NA)  # POINT warning

  # warning when geographic
  expect_warning(splm(y ~ x, exdata_sf_geo, spcov_type = spcov_type, estmethod = "reml"))
  expect_warning(splm(y ~ x, exdata_sf_NA, spcov_type = spcov_type, estmethod = "reml"), NA)
})

test_that("extra covr checks", {
  # random effects with cov_initial_search generics for 4 parameter families
  spcov_type <- "matern"
  # a random effect combined with matern's extra range parameter (and, below,
  # partitioning/anisotropy on top of that) leaves the automatic n <= 500
  # Satterthwaite ddf computation with a lot of covariance parameters to
  # estimate relative to this small fixture, which regularly makes its
  # covariance-parameter Hessian non-positive-definite -- spmodel already
  # warns and falls back to NULL ddf gracefully, so the warning is
  # expected/tolerated on the more heavily-parameterized calls below
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~group, partition_factor = ~group)), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~group, anisotropy = TRUE)), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group, partition_factor = ~group)), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group, anisotropy = TRUE)), NA)

  # random effects with cov_initial_search generics for 1 parameter families
  spcov_type <- "none"
  # same non-positive-definite-Hessian rationale as above
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group)), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group)), NA)

  # anisotropy for sv approaches
  spcov_initial_val <- spcov_initial(spcov_type = "exponential", rotate = 0.5, scale = 0.5, known = "given")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"), NA)

  # partition factors for sv approaches
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, spcov_type, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-wls", partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, spcov_type, xcoord = xcoord, ycoord = ycoord, estmethod = "sv-cl", partition_factor = ~group), NA)

  # more than one random effect
  expect_error(splm(y ~ x, exdata, "exponential", xcoord = xcoord, ycoord = ycoord, random = ~ group + subgroup), NA)
  # two random effects leave the automatic Satterthwaite ddf computation with
  # a lot of covariance parameters to estimate relative to this small
  # fixture, which regularly makes its covariance-parameter Hessian
  # non-positive-definite -- spmodel already warns and falls back to NULL
  # ddf gracefully, so the warning is expected/tolerated here
  expect_error(suppressWarnings(splm(y ~ x, exdata, "matern", xcoord = xcoord, ycoord = ycoord, random = ~ group + subgroup)), NA)
  expect_error(suppressWarnings(splm(y ~ x, exdata, "none", xcoord = xcoord, ycoord = ycoord, random = ~ group + subgroup)), NA)

  # anisotropy resets itself based on spcov initial
  spcov_initial_val <- spcov_initial(spcov_type = "exponential", rotate = 0, scale = 1, known = "given")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, anisotropy = TRUE), NA)

  # var adjust with anisotropy and random effects and partition factors
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml",
    anisotropy = TRUE, local = list(var_adjust = "theoretical")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml",
    random = ~group, local = list(var_adjust = "theoretical")
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml",
    partition_factor = ~group, local = list(var_adjust = "theoretical")
  ), NA)
})


test_that("examples run", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_output(print(spmod))
  expect_error(summary(spmod), NA)
  expect_error(tidy(spmod), NA)
  # different estimation method
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, estmethod = "sv-wls"), NA)
  # anisotropy
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, anisotropy = TRUE), NA)
  # random effects
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, random = ~group), NA)
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, random = ~ (x | group) + (x | subgroup)), NA)
  # partition factor
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, partition_factor = ~group), NA)
  # big data
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, local = TRUE), NA)
  ## parallel
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, local = list(parallel = TRUE, ncores = 2)), NA)
  # combining
  # optim non-convergence warning expected for this numerically marginal fixture
  expect_error(suppressWarnings(splm(y ~ x, exdata, "exponential", xcoord, ycoord, anisotropy = TRUE, random = ~group, partition_factor = ~group, local = TRUE)), NA)
  # spcov_initial
  spcov_initial_val <- spcov_initial("exponential", ie = 0, known = "ie")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord), NA)
  # randcov_initial ("group" is shorthand for "1 | group"))
  randcov_initial_val <- randcov_initial(1, nm = "group", known = "group")
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, random = ~group, randcov_initial = randcov_initial_val), NA)
})

test_that("errors occur", {
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, partition_factor = ~ group + subgroup))
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, local = list(method = "xyz")))
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, local = list(var_adjust = "xyz")))
  exdata2 <- exdata
  exdata2[1, "xcoord"] <- NA
  expect_error(splm(y ~ x, exdata2, "exponential", xcoord, ycoord))
  expect_error(suppressWarnings(splm(y ~ as.factor(x) + group, exdata, "exponential", xcoord, ycoord)))
  expect_error(splm(y ~ x, exdata, "exponential", xcoord = ycoord), NA) # changing to ycoord2 works
  expect_error(splm(y ~ x, exdata, "exponential", ycoord = xcoord))
  expect_error(splm(y ~ x, exdata, "xyz", xcoord, ycoord))
  # spcov_initial_val <- spcov_initial("xyz")
  # expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord))
  expect_error(splm(y ~ x, exdata, "car", xcoord, ycoord))
  spcov_initial_val <- spcov_initial("car")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord))
  expect_error(splm(y ~ x, exdata, "sar", xcoord, ycoord))
  spcov_initial_val <- spcov_initial("sar")
  expect_error(splm(y ~ x, exdata, spcov_initial = spcov_initial_val, xcoord = xcoord, ycoord = ycoord))
  exdata3 <- exdata
  exdata3$y <- as.character(exdata3$y)
  expect_error(splm(y ~ x, exdata3, "exponential", xcoord, ycoord))
  expect_error(splm(as.character(y) ~ x, exdata, "exponential", xcoord, ycoord))
  exdata3$y <- as.factor(exdata3$y)
  expect_error(splm(y ~ x, exdata3, "exponential", xcoord, ycoord))
  expect_error(splm(as.factor(y) ~ x, exdata, "exponential", xcoord, ycoord))
  exdata3$xcoord <- as.character(exdata3$xcoord)
  expect_error(splm(y ~ x, exdata3, "exponential", xcoord, ycoord))
  exdata3$ycoord <- as.character(exdata3$ycoord)
  expect_error(splm(y ~ x, exdata3, "exponential", xcoord, ycoord))
  exdata3$xcoord <- as.numeric(exdata3$xcoord)
  expect_error(splm(y ~ x, exdata3, "exponential", xcoord, ycoord))
  expect_error(splm(y ~ x, exdata, "exponential", xyz, ycoord))
  expect_error(splm(y ~ x, exdata, "exponential", "xyz", ycoord))
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, xyz))
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, "xyz"))
  exdata4 <- exdata
  exdata4$x2 <- exdata4$x
  expect_error(suppressWarnings(splm(y ~ x + x2, exdata4, "exponential", xcoord, ycoord)))
  exdata4$x[1] <- NA
  expect_error(splm(y ~ x, exdata4, "exponential", xcoord, ycoord))

  # anisotropy
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, anisotropy = TRUE, estmethod = "sv-wls"))
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, anisotropy = TRUE, estmethod = "sv-cl"))
  spcov_initial_val <- spcov_initial("exponential", rotate = 2)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"))
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"))
  spcov_initial_val <- spcov_initial("exponential", scale = 1)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"))
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"))
  spcov_initial_val <- spcov_initial("exponential", rotate = 2, scale = 1)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls"))
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl"))
  # random effects
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, estmethod = "sv-wls", random = ~group))
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, estmethod = "sv-cl", random = ~group))
})

test_that("messages occur", {
  expect_message(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord))
  spcov_initial_val <- spcov_initial("exponential")
  expect_message(splm(y ~ x, exdata, "exponential", xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val))
})

test_that("quoting arguments works", {
  spmod1 <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  spmod1$call <- NULL
  spmod2 <- splm(y ~ x, exdata, "exponential", "xcoord", "ycoord")
  spmod2$call <- NULL
  expect_equal(spmod1, spmod2)
})

test_that("no variance error works", {
  exdata$novar <- 1
  expect_error(splm(novar ~ x, exdata, "exponential", xcoord, ycoord))
})

test_that("offset works", {
  exdata$offset <- 2
  exdata$y2 <- exdata$y - exdata$offset
  spmod1 <- splm(y ~ x + offset(offset), exdata, "exponential", xcoord, ycoord)
  spmod2 <- splm(y2 ~ x, exdata, "exponential", xcoord, ycoord)
  expect_equal(fitted(spmod1), fitted(spmod2) + exdata$offset)
})

test_that("the model runs for partition and random effect group if there is an extra factor present", {
  spcov_type <- "exponential"
  exdata$group2 <- factor(exdata$group)
  levels(exdata$group2) <- c(levels(exdata$group2), ".new_group_level")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", partition_factor = ~group2), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group2), NA)
  spmod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group2)
  expect_equal(0, unname(fitted(spmod, type = "randcov")[["1 | group2"]]["group2.new_group_level"]))
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group2,
                    partition_factor = ~group2), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", partition_factor = ~group2, local = TRUE), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group2, local = TRUE), NA)
  spmod <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group2, local = TRUE)
  expect_equal(0, unname(fitted(spmod, type = "randcov")[["1 | group2"]]["group2.new_group_level"]))
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group2,
                    partition_factor = ~group, local = TRUE), NA)
})

test_that("emmeans works", {
  spcov_type <- "exponential"
  spmod <- splm(y ~ x * group, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(as.matrix(model.frame(delete.response(terms(spmod)), spmod$obdata)), as.matrix(emmeans::recover_data(spmod)))
  expect_error(emmeans::emmeans(spmod, ~ group, by = "x"), NA)
})

test_that("emmeans works missing", {
  spcov_type <- "exponential"
  spmod <- splm(y ~ x * group, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(as.matrix(model.frame(delete.response(terms(spmod)), spmod$obdata)), as.matrix(emmeans::recover_data(spmod)))
  expect_error(emmeans::emmeans(spmod, ~ group, by = "x"), NA)
})

test_that("covmatrix errors properly", {
  spcov_type <- "exponential"
  spmod <- splm(y ~ x * group, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
  expect_error(covmatrix(spmod, newdata = NULL))
  expect_error(covmatrix(spmod, cov_type = "xyz"), NA) # when newdata not specified cov_type silently ignored
  expect_error(covmatrix(spmod, newdata = spmod$newdata, cov_type = "xyz"))
})

test_that("range_constrain works", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", range_constrain = TRUE), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", range_constrain = TRUE), NA)
})

test_that("robust semivariogram works", {
  spcov_type <- "exponential"
  spmod1 <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls")
  spmod2 <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", robust = FALSE)
  spmod3 <- splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls", robust = TRUE)
  expect_true(identical(spmod1$esv, spmod2$esv))
  expect_true(!identical(spmod1$esv, spmod3$esv))
})


test_that("optim non-convergence warning fires for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  expect_warning(
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", control = list(maxit = 1)),
    "did not converge"
  )
  expect_warning(
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential"),
    NA
  )
})

test_that("reported ie matches covmatrix() when the estimated ie lands at the numerical floor", {
  # caribou + anisotropy is a real fixture whose estimated ie genuinely lands
  # below 1e-4 * de (the floor spcov_matrix.*() applies internally when
  # building Sigma); before floor_estimated_ie() was wired into use_gloglik*(),
  # the *reported* ie silently disagreed with the ie actually used to fit the
  # model
  spmod <- splm(z ~ water + tarp,
                data = caribou,
                spcov_type = "exponential", xcoord = x, ycoord = y,
                anisotropy = TRUE
  )
  spcov_coefs <- coef(spmod, type = "spcov")

  expect_equal(round(spcov_coefs[["ie"]], digits = 8), round(1e-4 * spcov_coefs[["de"]], digits = 8))

  # de + ie is the same-location (distance = 0) variance, i.e. the diagonal of
  # covmatrix() -- if the reported ie didn't match what was actually used to
  # build Sigma, this would fail
  expect_equal(spcov_coefs[["de"]] + spcov_coefs[["ie"]], unname(diag(covmatrix(spmod))[1]))
})


test_that("plot() works for esv()/eacf() when the cutoff exceeds the data's spatial extent", {
  # a cutoff beyond every observed pairwise distance leaves the outermost
  # bins empty (np = 0, gamma/acov = NA); the default ylim previously came
  # from max()/min() without na.rm = TRUE, so plot.window() failed with
  # "need finite 'ylim' values" -- esv() always hit this (gamma is never
  # negative, so its ylim branch always fires), while eacf() only hit it
  # when every non-NA acov value shared the same sign
  pdf(NULL)
  on.exit(dev.off())

  expect_error(plot(esv(sulfate ~ 1, sulfate, cutoff = 1e7)), NA)
  expect_error(plot(esv(sulfate ~ 1, sulfate, cutoff = 1e7, robust = TRUE)), NA)
  expect_error(plot(eacf(sulfate ~ 1, sulfate, cutoff = 1e7)), NA)

  # force eacf()'s all-positive-acov branch (same NA pattern) to confirm its
  # latent version of the same bug is fixed too
  e <- eacf(sulfate ~ 1, sulfate, cutoff = 1e7)
  e$acov <- abs(e$acov)
  expect_error(plot(e), NA)

  # sanity: an ordinary (non-empty-bin) cutoff still plots fine
  expect_error(plot(esv(sulfate ~ 1, sulfate)), NA)
  expect_error(plot(eacf(sulfate ~ 1, sulfate)), NA)
})

test_that("esv()/eacf() formula supports . as shorthand for all predictors, excluding coordinates", {
  # esv()/eacf() build their own model frame/matrix directly (they don't go
  # through get_data_object_splm()), so . support needs its own coordinate/
  # geometry exclusion -- see expand_formula_dot() in esv()/eacf()
  set.seed(14)
  n <- 30
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  y <- 1 + 2 * x1 - x2 + rnorm(n, sd = 0.3)
  d <- data.frame(y = y, x1 = x1, x2 = x2, xcoord = xcoord, ycoord = ycoord)

  esv_dot <- esv(y ~ ., data = d, xcoord = xcoord, ycoord = ycoord)
  esv_explicit <- esv(y ~ x1 + x2, data = d, xcoord = xcoord, ycoord = ycoord)
  expect_equal(esv_dot$gamma, esv_explicit$gamma)

  eacf_dot <- eacf(y ~ ., data = d, xcoord = xcoord, ycoord = ycoord)
  eacf_explicit <- eacf(y ~ x1 + x2, data = d, xcoord = xcoord, ycoord = ycoord)
  expect_equal(eacf_dot$acov, eacf_explicit$acov)

  # a reserved coordinate column can still be used explicitly alongside .
  expect_error(esv(y ~ xcoord + ., data = d, xcoord = xcoord, ycoord = ycoord), NA)

  # sf input: . must exclude the coordinates derived from geometry
  d_sf <- sf::st_as_sf(d, coords = c("xcoord", "ycoord"))
  esv_sf_dot <- esv(y ~ ., data = d_sf)
  expect_equal(esv_sf_dot$gamma, esv_explicit$gamma)
  eacf_sf_dot <- eacf(y ~ ., data = d_sf)
  expect_equal(eacf_sf_dot$acov, eacf_explicit$acov)

  # dist_matrix supplied directly (no coordinate columns to exclude)
  dm <- spdist(d, "xcoord", "ycoord")
  expect_error(esv(y ~ ., data = d, dist_matrix = dm), NA)
  expect_error(eacf(y ~ ., data = d, dist_matrix = dm), NA)
})

test_that("local = 'covariance' neighbor selection ranks by |covariance|, not raw covariance", {
  # spcov_types with negative covariance lobes (e.g. wave) need |covariance|
  # ranking so a strongly negatively-correlated neighbor isn't passed over in
  # favor of a weakly positively-correlated one -- verified directly here
  # since the ranking itself is inline in get_pred_splm() rather than a
  # separately-callable helper
  cov_vec <- c(-8, -7, 6, 5, 4, 3, 2, 1, 0.5, 0.1)
  size <- 3
  n <- length(cov_vec)
  old_idx <- order(cov_vec)[seq(n, n - size + 1)] # pre-fix (raw) ranking
  new_idx <- order(abs(cov_vec))[seq(n, n - size + 1)] # post-fix (|.|) ranking
  expect_false(setequal(old_idx, new_idx))
  expect_setequal(cov_vec[new_idx], c(-8, -7, 6))

  # regression check: predict()/decorrelate() with local = "covariance" still
  # run correctly end-to-end for a negative-lobe spcov_type (wave)
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
  spmod_wave <- splm(y ~ x, exdata, spcov_type = "wave", xcoord = xcoord, ycoord = ycoord)
  pred <- predict(spmod_wave, newdata = newexdata, local = list(method = "covariance", size = 10))
  expect_false(anyNA(pred))
  expect_true(all(is.finite(pred)))
  pred_block <- predict(spmod_wave, newdata = newexdata, block = TRUE, local = list(method = "covariance", size = 10))
  expect_true(is.finite(pred_block))
})

test_that("predict() errors informatively when newdata is missing coordinate columns", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod1 <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")

  # sanity: predict() still works normally when both coordinate columns are present
  expect_vector(predict(spmod1, newdata = newexdata))

  newexdata_noxcoord <- newexdata
  newexdata_noxcoord$xcoord <- NULL
  expect_error(predict(spmod1, newdata = newexdata_noxcoord), "coordinate column")
  expect_error(predict(spmod1, newdata = newexdata_noxcoord, block = TRUE), "coordinate column")

  newexdata_noycoord <- newexdata
  newexdata_noycoord$ycoord <- NULL
  expect_error(predict(spmod1, newdata = newexdata_noycoord), "coordinate column")

  newexdata_nocoord <- newexdata_noxcoord
  newexdata_nocoord$ycoord <- NULL
  expect_error(predict(spmod1, newdata = newexdata_nocoord), "coordinate column")

  # a 1D covariance (e.g. triangular) only truly needs xcoord in newdata --
  # ycoord is auto-filled with 0 and should not be flagged as missing
  spmod1d <- splm(y ~ x, exdata, spcov_type = "triangular", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  expect_vector(predict(spmod1d, newdata = newexdata_noycoord))
  expect_error(predict(spmod1d, newdata = newexdata_noxcoord), "coordinate column")

  # generic backstop: spdist_vectors() itself also catches a missing coordinate
  expect_error(
    spdist_vectors(data.frame(a = 1:3), data.frame(a = 4:5, b = 6:7), xcoord = "a", ycoord = "b", dim_coords = 2),
    "Coordinate column"
  )
})

test_that("predict() errors informatively when newdata has NA in a random-slope-only covariate", {
  set.seed(11)
  n <- 60
  x1 <- rnorm(n)
  x2 <- rnorm(n) # used only as a random slope, not a fixed effect
  xcoord <- runif(n)
  ycoord <- runif(n)
  grp <- factor(sample(letters[1:5], n, replace = TRUE))
  y <- 1 + 2 * x1 + rep(rnorm(5), length.out = n)[as.integer(grp)] * x2 + rnorm(n, sd = 0.3)
  d <- data.frame(y = y, x1 = x1, x2 = x2, grp = grp, xcoord = xcoord, ycoord = ycoord)

  spmod_slope <- splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~ (x2 | grp))

  newdata_good <- data.frame(
    x1 = c(0.5, -0.5), x2 = c(0.2, -0.3), grp = factor(c("a", "b"), levels = levels(grp)),
    xcoord = c(0.3, 0.6), ycoord = c(0.3, 0.6)
  )

  # sanity: works normally when x2 has no NA
  expect_vector(predict(spmod_slope, newdata_good))

  newdata_na <- newdata_good
  newdata_na$x2[1] <- NA
  expect_error(predict(spmod_slope, newdata_na), "Cannot have NA values in predictors.")
  expect_error(predict(spmod_slope, newdata_na, se.fit = TRUE), "Cannot have NA values in predictors.")
  expect_error(predict(spmod_slope, newdata_na, local = TRUE), "Cannot have NA values in predictors.")
})

test_that("predict() errors informatively when newdata has NA in a random intercept or partition factor grouping column", {
  set.seed(11)
  n <- 60
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  grp <- factor(sample(letters[1:5], n, replace = TRUE))
  y <- 1 + 2 * x1 + rep(rnorm(5), length.out = n)[as.integer(grp)] + rnorm(n, sd = 0.3)
  d <- data.frame(y = y, x1 = x1, grp = grp, xcoord = xcoord, ycoord = ycoord)

  newdata_good <- data.frame(
    x1 = c(0.5, -0.5), grp = factor(c("a", "b"), levels = levels(grp)),
    xcoord = c(0.3, 0.6), ycoord = c(0.3, 0.6)
  )
  newdata_na <- newdata_good
  newdata_na$grp[1] <- NA

  # NA must error, even though a genuinely new/unseen (non-NA) level is
  # handled gracefully (silently treated as having zero random effect/
  # partition factor contribution) -- these are different situations and
  # only the former should error
  newdata_newlevel <- newdata_good
  newdata_newlevel$grp <- factor(c("z", "b"), levels = c(levels(grp), "z"))

  spmod_int <- splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~ (1 | grp))
  expect_vector(predict(spmod_int, newdata_good))
  expect_vector(predict(spmod_int, newdata_newlevel))
  expect_error(predict(spmod_int, newdata_na), "Cannot have NA values in predictors.")
  expect_error(predict(spmod_int, newdata_na, se.fit = TRUE), "Cannot have NA values in predictors.")
  expect_error(predict(spmod_int, newdata_na, local = TRUE), "Cannot have NA values in predictors.")
  expect_error(predict(spmod_int, newdata_na, interval = "prediction"), "Cannot have NA values in predictors.")

  spmod_part <- splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, partition_factor = ~grp)
  expect_vector(predict(spmod_part, newdata_good))
  expect_vector(predict(spmod_part, newdata_newlevel))
  expect_error(predict(spmod_part, newdata_na), "Cannot have NA values in predictors.")
})

test_that("splm() errors informatively when a formula/random/partition_factor variable is not in data", {
  set.seed(12)
  n <- 30
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  grp <- factor(sample(letters[1:3], n, replace = TRUE))
  y <- 1 + 2 * x1 + rnorm(n, sd = 0.3)
  d <- data.frame(y = y, x1 = x1, grp = grp, xcoord = xcoord, ycoord = ycoord)

  # a same-named object in the calling environment (but not in data) should
  # not be silently picked up via ordinary formula scoping -- it should error
  not_a_col <- rnorm(n)
  not_a_group <- factor(sample(letters[1:3], n, replace = TRUE))

  expect_error(splm(y ~ not_a_col, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "not_a_col.*not found in data")
  expect_error(splm(not_a_col ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "not_a_col.*not found in data")
  expect_error(splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~not_a_group), "not_a_group.*not found in data")
  expect_error(splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, partition_factor = ~not_a_group), "not_a_group.*not found in data")

  # sanity: valid calls (including transformed predictors) still work
  expect_s3_class(splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "splm")
  expect_s3_class(splm(y ~ poly(x1, 2), data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "splm")
  expect_s3_class(splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~grp), "splm")
})

test_that("splm() formula supports . as shorthand for all predictors, excluding coordinates", {
  set.seed(13)
  n <- 30
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  y <- 1 + 2 * x1 - x2 + rnorm(n, sd = 0.3)
  d <- data.frame(y = y, x1 = x1, x2 = x2, xcoord = xcoord, ycoord = ycoord)

  mod_dot <- splm(y ~ ., data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  mod_explicit <- splm(y ~ x1 + x2, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  # . must expand to the non-coordinate predictors only
  expect_equal(names(coef(mod_dot)), names(coef(mod_explicit)))
  expect_equal(unname(coef(mod_dot)), unname(coef(mod_explicit)))
  expect_false(any(c("xcoord", "ycoord") %in% names(coef(mod_dot))))

  # model.frame()/model.matrix() must reuse the same (expanded) formula, not
  # re-expand . against obdata (which would pull xcoord/ycoord back in)
  expect_false(any(c("xcoord", "ycoord") %in% colnames(model.matrix(mod_dot))))
  expect_false("." %in% all.vars(mod_dot$formula))

  # the literal "." is preserved cosmetically in the printed call
  expect_true("." %in% all.vars(mod_dot$call$formula))

  # a reserved coordinate column can still be used explicitly alongside .
  mod_trend <- splm(y ~ xcoord + ., data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_true("xcoord" %in% names(coef(mod_trend)))
  expect_true(all(c("x1", "x2") %in% names(coef(mod_trend))))

  # sf input: . must exclude the coordinates derived from geometry, for both
  # POINT and (centroid-derived) other geometries
  d_sf <- sf::st_as_sf(d, coords = c("xcoord", "ycoord"))
  mod_sf_dot <- splm(y ~ ., data = d_sf, spcov_type = "exponential")
  expect_equal(unname(coef(mod_sf_dot)), unname(coef(mod_explicit)))

  # predict() on a . fit still works and stays consistent with newdata
  d_na <- d
  d_na$y[1:3] <- NA
  mod_dot_na <- splm(y ~ ., data = d_na, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_vector(predict(mod_dot_na))

  # random/partition_factor have no "everything else" meaning and must reject .
  expect_error(
    splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~.),
    "not supported in random"
  )
  expect_error(
    splm(y ~ x1, data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, partition_factor = ~.),
    "not supported in partition_factor"
  )
})

test_that("a user-supplied local$index must match the length of the non-missing response", {
  # local$index labels each row of obdata (the data actually used for
  # fitting) with a partition assignment; a row with a missing response is
  # excluded from obdata entirely (it becomes a prediction location instead),
  # so local$index must already be sized to the non-missing response vector,
  # not the original (possibly larger) data -- a mismatch is rejected
  # outright rather than guessed at, since split.data.frame()/split() don't
  # themselves validate that the grouping vector's length matches the data
  set.seed(12)
  n <- 40
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  y <- 1 + 2 * x1 + rnorm(n, sd = 0.3)
  d <- data.frame(y = y, x1 = x1, xcoord = xcoord, ycoord = ycoord)

  na_rows <- c(3, 10, 17, 25, 33)
  d_na <- d
  d_na$y[na_rows] <- NA
  n_obs <- n - length(na_rows)

  full_index <- sample(1:4, n, replace = TRUE) # wrong length: matches original n, not n_obs
  observed_index <- which(!is.na(d_na$y))
  correct_index <- full_index[observed_index] # right length: matches n_obs

  expect_error(
    get_data_object_splm(
      formula = y ~ x1, data = d_na, spcov_initial = spcov_initial("exponential"),
      xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml", anisotropy = FALSE,
      random = NULL, randcov_initial = NULL, partition_factor = NULL,
      local = list(index = full_index), range_constrain = FALSE
    ),
    "local\\$index must have the same length"
  )
  expect_error(
    splm(y ~ x1, data = d_na, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = list(index = full_index)),
    "local\\$index must have the same length"
  )

  # a correctly-sized index (matching the non-missing response vector) works
  data_object <- get_data_object_splm(
    formula = y ~ x1, data = d_na, spcov_initial = spcov_initial("exponential"),
    xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml", anisotropy = FALSE,
    random = NULL, randcov_initial = NULL, partition_factor = NULL,
    local = list(index = correct_index), range_constrain = FALSE
  )
  expect_length(data_object$local_index, n_obs)
  expect_equal(as.vector(data_object$local_index), as.vector(correct_index))

  expect_s3_class(
    splm(y ~ x1, data = d_na, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = list(index = correct_index)),
    "splm"
  )
})

