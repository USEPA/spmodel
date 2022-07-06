load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))
load(system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

test_that("the simulation runs for exponential", {
  spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for exponential anisotropy", {
  spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord, anisotropy = TRUE))
})

test_that("the simulation runs for exponential random effect", {
  spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  randcov_params_val <- randcov_params(group = 1)
  expect_vector(sprnorm(spcov_params_val,
    data = exdata, xcoord = xcoord, ycoord = ycoord,
    random = ~group, randcov_params = randcov_params_val
  ))
})

test_that("the simulation runs for exponential partitioning", {
  spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord, partition_factor = ~group))
})

test_that("the simulation runs for exponential (matrix)", {
  spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  sprnorm_matrix <- sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord, samples = 5)
  expect_equal(nrow(sprnorm_matrix), nrow(exdata))
  expect_equal(ncol(sprnorm_matrix), 5)
})

test_that("the simulation runs for car", {
  spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 1 / 2, extra = 0)
  expect_vector(sprnorm(spcov_params_val, data = exdata_poly))
})

test_that("the simulation runs for car unconnected", {
  spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 1 / 2, extra = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata_Upoly))
})

test_that("the simulation runs for none random effect", {
  spcov_params_val <- spcov_params("none", ie = 1)
  randcov_params_val <- randcov_params(group = 1)
  expect_vector(sprnorm(spcov_params_val,
    data = exdata, xcoord = xcoord, ycoord = ycoord,
    random = ~group, randcov_params = randcov_params_val
  ))
})

test_that("the simulation runs for car random effects", {
  spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 1 / 2, extra = 0)
  expect_vector(sprnorm(spcov_params_val, data = exdata_poly))
  randcov_params_val <- randcov_params(group = 1)
  expect_vector(sprnorm(spcov_params_val,
    data = exdata_poly, xcoord = xcoord, ycoord = ycoord,
    random = ~group, randcov_params = randcov_params_val
  ))
})

test_that("the simulation runs for car partition factor", {
  spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 1 / 2, extra = 0)
  expect_vector(sprnorm(spcov_params_val, data = exdata_poly, partition_factor = ~group))
})




test_that("the simulation runs for spherical", {
  spcov_params_val <- spcov_params("spherical", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for gaussian", {
  spcov_params_val <- spcov_params("gaussian", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for triangular", {
  spcov_params_val <- spcov_params("triangular", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord))
})

test_that("the simulation runs for circular", {
  spcov_params_val <- spcov_params("circular", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for cubic", {
  spcov_params_val <- spcov_params("cubic", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for penta", {
  spcov_params_val <- spcov_params("penta", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for cosine", {
  spcov_params_val <- spcov_params("cosine", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord))
})

test_that("the simulation runs for wave", {
  spcov_params_val <- spcov_params("wave", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for jbessel", {
  spcov_params_val <- spcov_params("jbessel", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for gravity", {
  spcov_params_val <- spcov_params("gravity", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for rquad", {
  spcov_params_val <- spcov_params("rquad", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for magnetic", {
  spcov_params_val <- spcov_params("magnetic", de = 1, ie = 1, range = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for none", {
  spcov_params_val <- spcov_params("none", ie = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata, mean = 4))
})

test_that("the simulation runs for matern", {
  spcov_params_val <- spcov_params("matern", de = 1, ie = 1, range = 1, extra = 1 / 2)
  expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for cauchy", {
  spcov_params_val <- spcov_params("cauchy", de = 1, ie = 1, range = 1, extra = 1 / 2)
  expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for pexponential", {
  spcov_params_val <- spcov_params("pexponential", de = 1, ie = 1, range = 1, extra = 1 / 2)
  expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
})

test_that("the simulation runs for sar", {
  spcov_params_val <- spcov_params("sar", de = 1, ie = 1, range = 1 / 2, extra = 0)
  expect_vector(sprnorm(spcov_params_val, data = exdata_poly))
})

test_that("the simulation runs for sar unconnected", {
  spcov_params_val <- spcov_params("sar", de = 1, ie = 1, range = 1 / 2, extra = 1)
  expect_vector(sprnorm(spcov_params_val, data = exdata_Upoly))
})
