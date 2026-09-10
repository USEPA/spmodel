skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)
skip_if_not_installed("ranger")

set.seed(1)

load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

test_that("the model runs", {
  spcov_type <- "car"
  expect_true(inherits(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type), "spautorRF"))
  expect_true(inherits(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "ml"), "spautorRF"))
})

test_that("the model runs", {
  spcov_type <- "car"
  num.tree <- 499
  spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
  expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val), NA)
  expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
})

test_that("the model list runs", {
  spcov_type <- c("car", "sar")
  expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type), NA)
  expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "ml"), NA)
})

test_that("the model list runs", {
  spcov_type <- c("car", "sar")
  num.tree <- 499
  spcov_initial_val <- lapply(spcov_type, function(x) spcov_initial(spcov_type = x, de = 1, ie = 1, range = 0.5, known = "de"))
  expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val), NA)
  expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
})

test_that("the model runs none ie", {
  spcov_type <- "none"
  expect_true(inherits(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type), "splmRF"))
  expect_true(inherits(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "ml"), "splmRF"))
  spcov_type <- "ie"
  expect_true(inherits(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type), "splmRF"))
  expect_true(inherits(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "ml"), "splmRF"))
})
