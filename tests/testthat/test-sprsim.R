load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))
load(system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)

test_that("splm-type sims work", {
  expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprnorm(spcov_params_val, mean = rnorm(NROW(exdata)), samples = 3, data = exdata, xcoord = xcoord, ycoord = ycoord), "matrix"))
})
