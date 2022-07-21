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
  randcov_initial_val <- randcov_initial(group = 1)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~ group + subgroup, randcov_initial = randcov_initial_val), NA)
})

test_that("the model runs for exponential (random nested subgroup)", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~ group / subgroup), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", random = ~ group / subgroup), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~ group / subgroup), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~ group / subgroup), NA)
  randcov_initial_val <- randcov_initial(group = 1)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", random = ~ group / subgroup, randcov_initial = randcov_initial_val), NA)
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
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "reml", anisotropy = TRUE
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "ml", anisotropy = TRUE
  ), NA)
})

test_that("the model runs for and random effects", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", anisotropy = TRUE, random = ~group
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", anisotropy = TRUE, random = ~group
  ), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "reml", anisotropy = TRUE, random = ~group
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "ml", anisotropy = TRUE, random = ~group
  ), NA)
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
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "reml", anisotropy = TRUE, partition_factor = ~group
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "ml", anisotropy = TRUE, partition_factor = ~group
  ), NA)
})

test_that("the model runs for and random effects and partitioning", {
  spcov_type <- "exponential"
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml", anisotropy = TRUE, random = ~group, partition_factor = ~group
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "ml", anisotropy = TRUE, random = ~group, partition_factor = ~group
  ), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "reml",
    anisotropy = TRUE, random = ~group, partition_factor = ~group
  ), NA)
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_initial = spcov_initial_val, estmethod = "ml",
    anisotropy = TRUE, random = ~group, partition_factor = ~group
  ), NA)
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
  expect_error(splm(y ~ x, exdata,
    xcoord = xcoord, ycoord = ycoord,
    spcov_type = spcov_type, estmethod = "reml",
    random = ~group, partition_factor = ~group, anisotropy = TRUE, local = TRUE
  ), NA)
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
  expect_warning(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"))
  expect_warning(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"))
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
  expect_warning(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"))
  expect_warning(splm(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"))
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
  # some positive definite issues -- need to figure out how to set zeros appropriately
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

  # if (!requireNamespace("sf", quietly = TRUE)) { # requireNamespace checks if sf is installed
  #   # dummy test so Skip test message are not printed
  #   expect_true(TRUE)
  # } else {
  # point data
  exdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"))

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
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_type <- "none"
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  spcov_type <- "matern"
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splm(y ~ x, exdata_poly, spcov_type = spcov_type, estmethod = "sv-cl"), NA)

  # exdata_sp <- sf::as_Spatial(exdata_sf) now saved in extdata

  # spcov_type <- "exponential"
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "reml"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "ml"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  # spcov_type <- "none"
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "reml"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "ml"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  # spcov_type <- "matern"
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "reml"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "ml"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  # expect_error(splm(y ~ x, exdata_sp, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  #
  # # polygon data
  # # exdata_poly_sp <- sf::as_Spatial(exdata_poly)
  # spcov_type <- "exponential"
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "reml"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "ml"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  # spcov_type <- "none"
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "reml"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "ml"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
  # spcov_type <- "matern"
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "reml"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "ml"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  # expect_error(splm(y ~ x, exdata_poly_sp, spcov_type = spcov_type, estmethod = "sv-cl"), NA)

  # }
})

test_that("extra covr checks", {
  # random effects with cov_initial_search generics for 4 parameter families
  spcov_type <- "matern"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~group, partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~group, anisotropy = TRUE), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group, partition_factor = ~group), NA)
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group, anisotropy = TRUE), NA)

  # random effects with cov_initial_search generics for 1 parameter families
  spcov_type <- "none"
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", random = ~group), NA)
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", random = ~group), NA)

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
  expect_error(splm(y ~ x, exdata, "matern", xcoord = xcoord, ycoord = ycoord, random = ~ group + subgroup), NA)
  expect_error(splm(y ~ x, exdata, "none", xcoord = xcoord, ycoord = ycoord, random = ~ group + subgroup), NA)

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
  expect_error(splm(y ~ x, exdata, "exponential", xcoord, ycoord, anisotropy = TRUE, random = ~group, partition_factor = ~group, local = TRUE), NA)
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
  expect_error(splm(y ~ as.factor(x) + group, exdata, "exponential", xcoord, ycoord))
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
