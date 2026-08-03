skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)
skip_if_not_installed("ranger")

set.seed(1)

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

test_that("the model runs for exponential", {
  spcov_type <- "exponential"
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
})

test_that("the model runs", {
  spcov_type <- "exponential"
  num.tree <- 499
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 1, known = "de")
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", num.trees = num.tree), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", num.trees = num.tree), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", num.trees = num.tree), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl", num.trees = num.tree), NA)
})

test_that("the model list runs", {
  spcov_type <- c("exponential", "matern")
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml"), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml"), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-wls"), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "sv-cl"), NA)
})

test_that("the model runs", {
  spcov_type <- c("exponential", "matern")
  num.tree <- 499
  spcov_initial_val <- lapply(spcov_type, function(x) spcov_initial(spcov_type = x, de = 1, ie = 1, range = 1, known = "de"))
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "reml", num.trees = num.tree), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "ml", num.trees = num.tree), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-wls", num.trees = num.tree), NA)
  expect_error(splmRF(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spcov_initial_val, estmethod = "sv-cl", num.trees = num.tree), NA)
})
