skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)
# SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
exdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"), crs = 5070)
load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
newexdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"), crs = 5070)
load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("the model runs for exponential", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val), NA)
})

test_that("the model runs for exponential (partition group)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, partition_factor = ~ group), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, partition_factor = ~ group), NA)
})

test_that("the model runs for exponential (random group)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group), NA)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val), NA)
})

test_that("the model runs for exponential (random and subgroup)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group + subgroup), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group + subgroup), NA)
  randcov_params_val <- randcov_params(group = 1, subgroup = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val), NA)
})

test_that("the model runs for exponential (random nested subgroup)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group / subgroup), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group / subgroup), NA)
  randcov_params_val <- randcov_params(group = 1, "group:subgroup" = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val), NA)
})

test_that("the model runs for exponential (random and partitioning)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group, partition_factor = ~ group), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group, partition_factor = ~ group), NA)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val, partition_factor = ~ group), NA)
})


test_that("the model runs for anisotropy", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, anisotropy = TRUE), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, rotate = 0.5, scale = 0.5)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val), NA)
})

test_that("the model runs for anisotropy (random effects)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, anisotropy = TRUE, random = ~ group), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, rotate = 0.5, scale = 0.5)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group), NA)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val), NA)
})

test_that("the model runs for anisotropy (partition factor)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, anisotropy = TRUE, partition_factor = ~ group), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, rotate = 0.5, scale = 0.5)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, partition_factor = ~ group), NA)
})

test_that("the model runs for anisotropy (random effects and partition factor)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, anisotropy = TRUE, random = ~ group, partition_factor = ~ group), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, rotate = 0.5, scale = 0.5)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group, partition_factor = ~ group), NA)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val, partition_factor = ~ group), NA)
})

test_that("the model runs for missing data", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val), NA)
})

test_that("the model runs for missing data (random group)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group), NA)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val), NA)
})

test_that("the model runs for big data", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, local = TRUE), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, local = list(parallel = TRUE, ncores = 2)), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, local = list(size = 30)), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, local = list(size = 10)), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, local = list(size = 10, method = "distance")), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, local = list(size = 10, method = "all")), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group, local = TRUE), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, anisotropy = TRUE, local = TRUE), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, partition_factor = ~ group, local = TRUE), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val, local = TRUE), NA)
})

test_that("the model runs for different ordering", {
  spcov_type <- "exponential"
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, ordering = "grts"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, ordering = "grts"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, ordering = "random"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, ordering = "random"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, ordering = "none"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, ordering = "none"), NA)
})

test_that("the model runs for evaluate test", {
  spcov_type <- "exponential"
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = "x", ycoord = y, spcov_params = spcov_params_val, evaluate_test = TRUE), NA)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = "ycoord", spcov_params = spcov_params_val, randcov_params = randcov_params_val, evaluate_test = TRUE), NA)
})

test_that("the model runs for different coordinates", {
  spcov_type <- "exponential"
  exdata$xc <- exdata$xcoord
  exdata$yc <- exdata$ycoord
  newexdata$xc <- newexdata$xcoord
  newexdata$yc <- newexdata$ycoord
  expect_error(decorrelate(y ~ x, exdata, xcoord = xc, ycoord = "yc", spcov_type = spcov_type), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = "xc", ycoord = yc, spcov_type = spcov_type), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = "xc", ycoord = yc, spcov_params = spcov_params_val), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xc, ycoord = "yc", spcov_params = spcov_params_val), NA)
  mod <- decorrelate(y ~ x, exdata, xcoord = "xc", ycoord = "yc", spcov_type = spcov_type)
  expect_vector(predict(mod, newdata = newexdata))
})

test_that("the model runs for other spatial covariances", {
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "gaussian"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "triangular"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "circular"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "cubic"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "pentaspherical"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "cosine"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "wave"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "jbessel"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "gravity"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "rquad"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "magnetic"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "cauchy"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "pexponential"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none"), NA)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "ie"), NA)
})

test_that("the model runs for sf objects", {
  # point data
  # exdata_sf_geo <- sf::st_transform(exdata_sf, crs = 4326)
  # exdata_sf_NA <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"), crs = NA)
  expect_error(decorrelate(y ~ x, exdata_sf, spcov_type = "exponential"), NA)
  expect_warning(expect_warning(decorrelate(y ~ x, exdata_poly, spcov_type = "exponential")))
  preds <- predict(decorrelate(y ~ x, exdata_sf, spcov_type = "exponential"), newdata = newexdata_sf)
  expect_vector(preds)
})

# predict
test_that("prediction works", {

  spcov_type <- "exponential"

  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type)
  expect_error(predict(mod, newdata = newexdata), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val)
  expect_error(predict(mod, newdata = newexdata), NA)

  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group, partition_factor = ~ group, anisotropy = TRUE)
  expect_error(predict(mod, newdata = newexdata), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, rotate = 0.5, scale = 0.5)
  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group, partition_factor = ~ group)
  expect_error(predict(mod, newdata = newexdata), NA)
  randcov_params_val <- randcov_params(group = 1)
  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val, partition_factor = ~ group)
  expect_error(predict(mod, newdata = newexdata), NA)

  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, local = TRUE)
  expect_error(predict(mod, newdata = newexdata), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, local = TRUE)
  expect_error(predict(mod, newdata = newexdata), NA)

  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type)
  expect_error(predict(mod, newdata = newexdata, local = TRUE), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val)
  expect_error(predict(mod, newdata = newexdata, local = list(size = 20, method = "distance", parallel = TRUE, ncores = 2)), NA)

  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, local = TRUE)
  expect_error(predict(mod, newdata = newexdata, local = list(method = "all")), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  mod <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, local = TRUE)
  expect_error(predict(mod, newdata = newexdata, local = TRUE), NA)

})

test_that("direct functions work", {
  spcov_type <- "exponential"
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  mod <- decorrelate_data(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val)
  expect_s3_class(mod, "decorrelate_data")
  newmod <- decorrelate_newdata(mod, newdata = newexdata)
  expect_s3_class(newmod, "decorrelate_newdata")
  preds <- rnorm(NROW(newexdata))
  newpreds <- recorrelate_newdata(newmod, preds)
  expect_vector(newpreds)

  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, rotate = 0.5, scale = 0.5)
  randcov_params_val <- c(group = 1)
  mod <- decorrelate_data(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val, partition_factor = ~ group, local = TRUE, ordering = "random")
  expect_s3_class(mod, "decorrelate_data")
  newmod <- decorrelate_newdata(mod, newdata = newexdata, local = TRUE)
  expect_s3_class(newmod, "decorrelate_newdata")
  preds <- rnorm(NROW(newexdata))
  newpreds <- recorrelate_newdata(newmod, preds)
  expect_vector(newpreds)

  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  mod <- decorrelate_data(y ~ x, exdata_sf, spcov_params = spcov_params_val)
  expect_s3_class(mod, "decorrelate_data")
  newmod <- decorrelate_newdata(mod, newdata = newexdata_sf)
  expect_s3_class(newmod, "decorrelate_newdata")
  preds <- rnorm(NROW(newexdata_sf))
  newpreds <- recorrelate_newdata(newmod, preds)
  expect_vector(newpreds)

})

test_that("direct functions work (different coordinates)", {
  spcov_type <- "exponential"
  exdata$xc <- exdata$xcoord
  exdata$yc <- exdata$ycoord
  newexdata$xc <- newexdata$xcoord
  newexdata$yc <- newexdata$ycoord
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  mod <- decorrelate_data(y ~ x, exdata, xcoord = "xc", ycoord = yc, spcov_params = spcov_params_val)
  expect_s3_class(mod, "decorrelate_data")
  newmod <- decorrelate_newdata(mod, newdata = newexdata)
  expect_s3_class(newmod, "decorrelate_newdata")
  preds <- rnorm(NROW(newexdata))
  newpreds <- recorrelate_newdata(newmod, preds)
  expect_vector(newpreds)
})

test_that("grid works", {

  spcov_type <- "exponential"
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  randcov_params_val <- c(group = 1)

  expect_error(decorrelate_grid(y ~ x, exdata, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord), NA)
  expect_error(decorrelate_grid(y ~ x, exdata, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, anisotropy = TRUE), NA)
  expect_error(decorrelate_grid(y ~ x, exdata, spcov_params = spcov_params_val, xcoord = xcoord, ycoord = ycoord), NA)

  expect_error(decorrelate_grid(y ~ x, exdata, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, random = ~ group), NA)
  expect_error(decorrelate_grid(y ~ x, exdata, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, random = ~ group, anisotropy = TRUE), NA)
  expect_error(decorrelate_grid(y ~ x, exdata, spcov_params = spcov_params_val, xcoord = xcoord, ycoord = ycoord, random = ~ group), NA)
  expect_error(decorrelate_grid(y ~ x, exdata, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, randcov_params = randcov_params_val), NA)
  expect_error(decorrelate_grid(y ~ x, exdata, spcov_params = spcov_params_val, xcoord = xcoord, ycoord = ycoord, randcov_params = randcov_params_val), NA)

  expect_error(decorrelate_grid(y ~ x, exdata_sf, spcov_type = spcov_type), NA)
  expect_error(decorrelate_grid(y ~ x, exdata_sf, spcov_type = spcov_type, anisotropy = TRUE), NA)
  expect_error(decorrelate_grid(y ~ x, exdata_sf, spcov_type = spcov_type, random = ~ group), NA)
  expect_error(decorrelate_grid(y ~ x, exdata_sf, spcov_params = spcov_params_val), NA)

})

test_that("grid works (different coordinates)", {
  exdata$xc <- exdata$xcoord
  exdata$yc <- exdata$ycoord
  spcov_type <- "exponential"
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)

  expect_error(decorrelate_grid(y ~ x, exdata, spcov_type = spcov_type, xcoord = xc, ycoord = yc), NA)
  expect_error(decorrelate_grid(y ~ x, exdata, spcov_type = spcov_type, xc = xcoord, ycoord = "yc", anisotropy = TRUE), NA)
  expect_error(decorrelate_grid(y ~ x, exdata, spcov_params = spcov_params_val, xcoord = "xc", ycoord = "yc"), NA)

})


# training list methods
test_that("the model runs for exponential", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, training = list(method = "cv")), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, training = list(method = "cv")), NA)
})

test_that("the model runs for exponential (partition group)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, partition_factor = ~ group, training = list(method = "cv")), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, partition_factor = ~ group, training = list(method = "cv")), NA)
})

test_that("the model runs for exponential (random group)", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group, training = list(method = "cv")), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group, training = list(method = "cv")), NA)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val, training = list(method = "cv")), NA)
})

test_that("the model runs for matern and local", {
  spcov_type <- "matern"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, training = list(method = "cv"), local = TRUE), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, extra = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, training = list(method = "cv"), local = TRUE), NA)
})

test_that("the model runs for matern and local (partition group)", {
  spcov_type <- "matern"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, partition_factor = ~ group, training = list(method = "cv"), local = TRUE), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, extra = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, partition_factor = ~ group, training = list(method = "cv"), local = TRUE), NA)
})

test_that("the model runs for matern and local (random group)", {
  spcov_type <- "matern"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group, training = list(method = "cv"), local = TRUE), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1, extra = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group, training = list(method = "cv"), local = TRUE), NA)
  randcov_params_val <- randcov_params(group = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val, training = list(method = "cv"), local = TRUE), NA)
})

# dense grid methods
test_that("the model runs for exponential", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, dense_grid = FALSE), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val,  dense_grid = FALSE, training = list(method = "cv")), NA)
})

test_that("the model runs for variables and local", {
  spcov_type <- "exponential"
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, random = ~ group, partition_factor = ~ group, anisotropy = TRUE, dense_grid = FALSE), NA)
  spcov_params_val <- spcov_params(spcov_type = spcov_type, de = 1, ie = 1, range = 1)
  expect_error(decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, random = ~ group, partition_factor = ~ group, dense_grid = FALSE, training = list(method = "cv")), NA)
})

test_that("decorrelate() formula supports . as shorthand for all predictors, excluding coordinates", {
  # decorrelate()/decorrelate_data() build their data object via the same
  # get_data_object_splm() used by splm(), so . support (and its exclusion of
  # xcoord/ycoord) is inherited automatically -- this just confirms that
  # end-to-end, including through decorrelate_newdata()'s reuse of the
  # already-resolved terms object
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
  # exdata also has group/subgroup columns (used elsewhere for random effect
  # tests); restrict to x/xcoord/ycoord so . unambiguously matches ~ x below
  exdata <- exdata[, c("y", "x", "xcoord", "ycoord")]

  params <- spcov_params("exponential", de = 1, ie = 0.2, range = 1e5)

  decorr_dot <- decorrelate_data(y ~ ., data = exdata, xcoord = xcoord, ycoord = "ycoord", spcov_params = params)
  decorr_explicit <- decorrelate_data(y ~ x, data = exdata, xcoord = xcoord, ycoord = "ycoord", spcov_params = params)
  expect_equal(colnames(decorr_dot$X), colnames(decorr_explicit$X))
  expect_equal(decorr_dot$X, decorr_explicit$X)
  expect_equal(decorr_dot$tX, decorr_explicit$tX)
  expect_false(any(c("xcoord", "ycoord") %in% colnames(decorr_dot$X)))

  decorr_full <- decorrelate(y ~ ., exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = "ycoord")
  preds_dot <- predict(decorr_full, newdata = newexdata)
  expect_vector(preds_dot)
})

test_that("decorrelate() does not produce NaN with a nested random effect", {
  # a nested random effect (random = ~ group / subgroup) makes a neighbor
  # near-perfectly predictive of some observations, which can push the
  # conditional variance w = 1 - r0'SigmaInv r0 just below 0 due to floating
  # point roundoff -- sqrt(w) then produced NaN before w was floored at 0
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  expect_warning(
    decorr1 <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", random = ~ group / subgroup),
    NA
  )
  expect_false(anyNA(decorr1$decorrelate_data$tX))
  expect_false(anyNA(decorr1$decorrelate_data$ty))
  expect_true(all(is.finite(as.matrix(decorr1$decorrelate_data$tX))))
  expect_true(all(is.finite(decorr1$decorrelate_data$ty)))

  spcov_params_val <- spcov_params(spcov_type = "exponential", de = 1, ie = 1, range = 1)
  randcov_params_val <- randcov_params(group = 1, "group:subgroup" = 1)
  expect_warning(
    decorr2 <- decorrelate(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_params = spcov_params_val, randcov_params = randcov_params_val),
    NA
  )

  exdata_na <- exdata
  exdata_na$y[1:5] <- NA
  expect_warning(
    decorr3 <- decorrelate(y ~ x, exdata_na, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", random = ~ group / subgroup),
    NA
  )
  preds3 <- predict(decorr3)
  expect_false(anyNA(preds3))
  expect_true(all(is.finite(preds3)))
})