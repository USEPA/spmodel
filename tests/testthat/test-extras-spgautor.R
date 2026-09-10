skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

# SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

# add variables
n <- NROW(exdata_poly)
exdata_poly$bern <- rbinom(n, size = 1, prob = 0.5)
exdata_poly$bernfac <- factor(ifelse(exdata_poly$bern == 0, "a", "b"))
exdata_poly$size <- 5
exdata_poly$bin <- rbinom(n, size = exdata_poly$size, prob = 0.5)
exdata_poly$prop <- runif(n)
exdata_poly$count <- rpois(n, lambda = 5)
exdata_poly$cont <- rgamma(n, shape = 1, rate = 1)
exdata_poly$offset <- 2

# save W and M
W <- sf::st_intersects(exdata_poly, sparse = FALSE)
diag(W) <- 0
W <- 1 * Matrix::Matrix(W, sparse = TRUE)
W_rowsums <- Matrix::rowSums(W)
Wbar <- W / W_rowsums
M <- Matrix(diag(NROW(W))) # for row-unstandardized / as matrix
Mvec <- 1 / W_rowsums # for row-standardized / as vector

test_that("the model runs for binomial data", {
  spgmod <- spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml")
  expect_s3_class(spgmod, "spgautor")
  expect_vector(AUROC(spgmod))
  expect_error(spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml", range_positive = FALSE), NA)
  expect_error(spgautor(bernfac ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "ml"), NA)
  # perfect separation warning expected for this simulated binomial fixture
  expect_error(suppressWarnings(spgautor(cbind(bin, size) ~ x, family = "binomial", data = exdata_poly, spcov_type = "sar", estmethod = "reml")), NA)
  # causes an error with ml estimation as de is near zero and ie is zero, which makes inverse unstable
  # need to implement a diagonal tolerance for gautor models
  # expect_error(spgautor(cbind(bin, size) ~ x, family = "binomial", data = exdata_poly, spcov_type = "sar", estmethod = "reml"), NA)
  spgmod <- spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml")
  expect_vector(AUROC(spgmod))

  # complicated models
  expect_error(spgautor(bern ~ x,
    family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml",
    W = W, row_st = FALSE, M = M
  ), NA)
  expect_error(spgautor(bernfac ~ x,
    family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml",
    W = W, row_st = FALSE, M = M
  ), NA)
  # perfect separation warning expected for this simulated binomial fixture
  expect_error(suppressWarnings(spgautor(cbind(bin, size) ~ x + offset(offset),
    family = "binomial", data = exdata_poly, spcov_type = "sar", estmethod = "ml",
    random = ~group
  )), NA)
  expect_error(spgautor(y > 0 ~ x,
    family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml",
    W = W, row_st = FALSE, M = M
  ), NA)
})

test_that("the model runs for proportion data", {
  expect_error(spgautor(prop ~ x, family = "beta", data = exdata_poly, spcov_type = "car", estmethod = "reml"), NA)
  expect_error(spgautor(prop ~ x, family = "beta", data = exdata_poly, spcov_type = "car", estmethod = "reml", range_positive = FALSE), NA)
  expect_error(spgautor(prop ~ x, family = beta, data = exdata_poly, spcov_type = "sar", estmethod = "reml"), NA)

  # complicated models
  expect_error(spgautor(prop ~ x,
    family = "beta", data = exdata_poly, spcov_type = "car", estmethod = "reml",
    random = ~ group + subgroup
  ), NA)
  expect_error(spgautor(prop ~ x + offset(offset),
    family = beta, data = exdata_poly, spcov_type = "sar", estmethod = "reml",
    W = W, partition_factor = ~group
  ), NA)
})

test_that("the model runs for count data", {
  spgmod <- spgautor(count ~ x, family = poisson, data = exdata_poly, spcov_type = "sar", estmethod = "reml")
  expect_s3_class(spgmod, "spgautor")
  expect_error(AUROC(spgmod))
  expect_error(spgautor(count ~ x, family = "nbinomial", data = exdata_poly, spcov_type = "car", estmethod = "ml"), NA)
  expect_error(spgautor(count ~ x, family = "nbinomial", data = exdata_poly, spcov_type = "car", estmethod = "ml", range_positive = FALSE), NA)

  # complicated models
  expect_error(spgautor(count ~ x,
    family = poisson, data = exdata_poly, spcov_type = "sar", estmethod = "reml",
    random = ~ (x | subgroup)
  ), NA)
  expect_error(spgautor(count ~ x + offset(offset),
    family = "nbinomial", data = exdata_poly, spcov_type = "car", estmethod = "ml",
    W = Wbar, row_st = FALSE, M = Mvec
  ), NA)
})

test_that("the model runs for continuous data", {
  expect_error(spgautor(cont ~ x, family = "Gamma", data = exdata_poly, spcov_type = "car", estmethod = "reml"), NA)
  expect_error(spgautor(cont ~ x, family = inverse.gaussian, data = exdata_poly, spcov_type = "sar", estmethod = "ml"), NA)
  expect_error(spgautor(cont ~ x, family = inverse.gaussian, data = exdata_poly, spcov_type = "car", estmethod = "reml", range_positive = FALSE), NA)
  # SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
  expect_error(spgautor(cont ~ x, family = gaussian, data = exdata_poly, spcov_type = "car", estmethod = "reml"))

  # complicated models
  expect_error(spgautor(cont ~ x + offset(offset),
    family = "Gamma", data = exdata_poly, spcov_type = "car", estmethod = "reml",
    W = Wbar, row_st = FALSE, M = Mvec
  ), NA)
  expect_error(spgautor(cont ~ x,
    family = inverse.gaussian, data = exdata_poly, spcov_type = "sar", estmethod = "ml",
    random = ~subgroup, partition_factor = ~group
  ), NA)
  ## SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
  expect_error(spgautor(cont ~ x,
    family = gaussian, data = exdata_poly, spcov_type = "car", estmethod = "reml",
    random = ~subgroup
  ))
})

test_that("the model runs on other data sets", {
  expect_error(spgautor(abs(y) ~ x, family = "Gamma", data = exdata_Mpoly, spcov_type = "car", estmethod = "reml"), NA)
  expect_error(spgautor(abs(y) ~ x, family = Gamma, data = exdata_Upoly, spcov_type = "sar", estmethod = "reml"), NA)

  # complicated models
  expect_error(spgautor(abs(y) ~ x,
    family = "Gamma", data = exdata_Mpoly, spcov_type = "car", estmethod = "reml",
    W = W, partition_factor = ~subgroup
  ), NA)

  W <- sf::st_intersects(exdata_Upoly, sparse = FALSE)
  diag(W) <- 0
  W <- 1 * Matrix::Matrix(W, sparse = TRUE)
  expect_error(spgautor(abs(y) ~ x,
    family = Gamma, data = exdata_Upoly, spcov_type = "sar", estmethod = "reml",
    W = W, random = ~subgroup
  ), NA)
})

test_that("emmeans works", {
  spcov_type <- "car"
  spgmod <- spgautor(abs(y) ~ x * group, family = "Gamma", exdata_poly, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(as.matrix(model.frame(delete.response(terms(spgmod)), spgmod$data[spgmod$observed_index, , drop = FALSE])), as.matrix(emmeans::recover_data(spgmod)))
  expect_error(emmeans::emmeans(spgmod, ~ group, by = "x"), NA)
})

test_that("emmeans works missing", {
  spcov_type <- "car"
  spgmod <- spgautor(abs(y) ~ x * group, family = "Gamma", exdata_Mpoly, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(as.matrix(model.frame(delete.response(terms(spgmod)), spgmod$data[spgmod$observed_index, , drop = FALSE])), as.matrix(emmeans::recover_data(spgmod)))
  expect_error(emmeans::emmeans(spgmod, ~ group, by = "x"), NA)
})

test_that("point distance works missing", {
  exdata_sf <- st_as_sf(exdata, coords = c("xcoord", "ycoord"), crs = NA)
  spcov_type <- "sar"
  expect_error(spgautor(abs(y) ~ x, family = "Gamma", exdata_sf, spcov_type = spcov_type, cutoff = 1), NA)
  expect_error(spgautor(abs(y) ~ x, family = Gamma, exdata_sf, spcov_type = spcov_type, cutoff = 1, row_st = FALSE), NA)
  expect_error(spgautor(abs(y) ~ x, family = "Gamma", exdata_sf, spcov_type = spcov_type, cutoff = NULL)) # can't be NULL
  expect_error(spgautor(abs(y) ~ x, family = Gamma, exdata_sf, spcov_type = spcov_type, cutoff = 1e-8)) # too small of distance so no neighbors
})

test_that("spcov_type none and ie work properly", {
  spcov_type <- "none"
  spgmod <- spgautor(bern ~ x * group, family = "binomial", exdata_poly, spcov_type = spcov_type, estmethod = "reml")
  expect_true(inherits(spgmod, "spglm"))
  spgmod_lm <- spglm(bern ~ x * group, family = "binomial", exdata_poly, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(spgmod, spgmod_lm)
  spcov_type <- "ie"
  spgmod <- spgautor(bern ~ x * group, family = binomial, exdata_poly, spcov_type = spcov_type, estmethod = "reml")
  expect_true(inherits(spgmod, "spglm"))
  spgmod_lm <- spglm(bern ~ x * group, family = binomial, exdata_poly, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(spgmod, spgmod_lm)
  spgmod <- spgautor(bern ~ x * group, family = "binomial", exdata_poly, spcov_type = c("none", "ie", "car", "sar"), estmethod = "reml")
  expect_true(inherits(spgmod, "spgautor_list"))
  expect_true(inherits(spgmod$none, "spglm"))
  expect_true(inherits(coef(spgmod$none, type = "spcov"), "none"))
  expect_true(inherits(spgmod$ie, "spglm"))
  expect_true(inherits(coef(spgmod$ie, type = "spcov"), "ie"))
  expect_true(inherits(spgmod$car, "spgautor"))
  expect_true(inherits(coef(spgmod$car, type = "spcov"), "car"))
  expect_true(inherits(spgmod$sar, "spgautor"))
  expect_true(inherits(coef(spgmod$sar, type = "spcov"), "sar"))
  expect_error(glances(spgmod), NA)
})

test_that("optim non-convergence warning fires for spgautor", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  expect_warning(
    spgautor(abs(y) ~ x, family = Gamma, exdata_poly, spcov_type = "car", estmethod = "reml", control = list(maxit = 1)),
    "did not converge"
  )
  expect_warning(
    spgautor(abs(y) ~ x, family = Gamma, exdata_poly, spcov_type = "car", estmethod = "reml"),
    NA
  )
})

test_that("spgautor() errors informatively when a formula/random/partition_factor variable is not in data", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  # a same-named object in the calling environment (but not in data) should
  # not be silently picked up via ordinary formula scoping -- it should error
  not_a_col <- rnorm(NROW(exdata_poly))
  not_a_group <- factor(sample(letters[1:3], NROW(exdata_poly), replace = TRUE))

  expect_error(spgautor(abs(y) ~ not_a_col, family = "Gamma", data = exdata_poly, spcov_type = "car"), "not_a_col.*not found in data")
  expect_error(spgautor(not_a_col ~ x, family = "Gamma", data = exdata_poly, spcov_type = "car"), "not_a_col.*not found in data")
  expect_error(spgautor(abs(y) ~ x, family = "Gamma", data = exdata_poly, spcov_type = "car", random = ~not_a_group), "not_a_group.*not found in data")
  expect_error(spgautor(abs(y) ~ x, family = "Gamma", data = exdata_poly, spcov_type = "car", partition_factor = ~not_a_group), "not_a_group.*not found in data")

  # sanity: a valid call still works
  expect_s3_class(spgautor(abs(y) ~ x, family = "Gamma", data = exdata_poly, spcov_type = "car"), "spgautor")
})

test_that("spgautor() formula supports . as shorthand for all predictors, excluding the sf geometry column", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  d <- exdata_poly[, c("y", "x")] # sf's `[` keeps geometry regardless of selection

  mod_dot <- spgautor(abs(y) ~ ., family = "Gamma", data = d, spcov_type = "car")
  mod_explicit <- spgautor(abs(y) ~ x, family = "Gamma", data = d, spcov_type = "car")
  expect_equal(names(coef(mod_dot)), names(coef(mod_explicit)))
  expect_equal(unname(coef(mod_dot)), unname(coef(mod_explicit)))
  expect_false("geometry" %in% colnames(model.matrix(mod_dot)))

  expect_error(
    spgautor(abs(y) ~ x, family = "Gamma", data = d, spcov_type = "car", random = ~.),
    "not supported in random"
  )
})

test_that("predict() only allows newdata = object$newdata for spgautor()", {
  # spgautor() prediction locations are fixed when the model is fit (they
  # determine the neighbor structure W/M used throughout fitting), so a
  # different newdata cannot be honored at predict() time
  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  gamod <- spgautor(abs(y) ~ x, family = "Gamma", exdata_Mpoly, spcov_type = "car")

  expect_vector(predict(gamod))
  expect_equal(predict(gamod), predict(gamod, newdata = gamod$newdata))

  modified_newdata <- gamod$newdata
  modified_newdata$x <- modified_newdata$x + 1
  expect_error(predict(gamod, newdata = modified_newdata), "newdata cannot be specified")
  expect_error(predict(gamod, newdata = exdata_poly), "newdata cannot be specified")

  # a model with no missing data at all should still error informatively
  gamod_full <- spgautor(abs(y) ~ x, family = "Gamma", exdata_poly, spcov_type = "car")
  expect_error(predict(gamod_full), "No missing data to predict")
})