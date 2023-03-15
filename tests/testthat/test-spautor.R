set.seed(1)

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

#### CRAN check
test_that("the model runs for connected sites", {
  expect_error(spautor(y ~ x, exdata_poly, spcov_type = "car"), NA)
})

test_that("the model runs for unconnected sites", {
  expect_error(spautor(y ~ x, exdata_Upoly, spcov_type = "sar"), NA)
})
