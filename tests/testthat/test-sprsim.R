test_that("splm and spglm type sims work", {

  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  load(system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

  spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)

  # sprnorm
  expect_vector(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprnorm(spcov_params_val, mean = rnorm(NROW(exdata)), samples = 3, data = exdata, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprpois
  expect_vector(sprpois(spcov_params_val, data = exdata_poly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprpois(spcov_params_val, mean = rnorm(NROW(exdata)), samples = 3, data = exdata, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprnbinom
  expect_vector(sprnbinom(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprnbinom(spcov_params_val, mean = rnorm(NROW(exdata)), samples = 3, data = exdata, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprbinom
  expect_vector(sprbinom(spcov_params_val, data = exdata_Upoly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprbinom(spcov_params_val, mean = rnorm(NROW(exdata)), samples = 3, data = exdata, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprbeta
  expect_vector(sprbeta(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprbeta(spcov_params_val, mean = rnorm(NROW(exdata)), samples = 3, data = exdata, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprgamma
  expect_vector(sprgamma(spcov_params_val, data = exdata_poly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprgamma(spcov_params_val, mean = rnorm(NROW(exdata)), samples = 3, data = exdata, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprinvgauss
  expect_vector(sprinvgauss(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprinvgauss(spcov_params_val, mean = rnorm(NROW(exdata_Upoly)), samples = 3, data = exdata_Upoly, xcoord = xcoord, ycoord = ycoord), "matrix"))
})

test_that("spautor and spgautor type sims work", {

  load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  load(system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

  spcov_params_val <- spcov_params("car", de = 1, ie = 0, extra = 1, range = 0.2)

  # sprnorm
  expect_vector(sprnorm(spcov_params_val, data = exdata_poly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprnorm(spcov_params_val, mean = rnorm(NROW(exdata_poly)), samples = 3, data = exdata_poly, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprpois
  expect_vector(sprpois(spcov_params_val, data = exdata_poly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprpois(spcov_params_val, mean = rnorm(NROW(exdata_poly)), samples = 3, data = exdata_poly, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprnbinom
  expect_vector(sprnbinom(spcov_params_val, data = exdata_poly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprnbinom(spcov_params_val, mean = rnorm(NROW(exdata_poly)), samples = 3, data = exdata_poly, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprbinom
  expect_vector(sprbinom(spcov_params_val, data = exdata_Upoly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprbinom(spcov_params_val, mean = rnorm(NROW(exdata_poly)), samples = 3, data = exdata_poly, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprbeta
  expect_vector(sprbeta(spcov_params_val, data = exdata_poly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprbeta(spcov_params_val, mean = rnorm(NROW(exdata_poly)), samples = 3, data = exdata_poly, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprgamma
  expect_vector(sprgamma(spcov_params_val, data = exdata_poly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprgamma(spcov_params_val, mean = rnorm(NROW(exdata_poly)), samples = 3, data = exdata_poly, xcoord = xcoord, ycoord = ycoord), "matrix"))

  # sprinvgauss
  expect_vector(sprinvgauss(spcov_params_val, data = exdata_poly, xcoord = xcoord, ycoord = ycoord))
  expect_true(inherits(sprinvgauss(spcov_params_val, mean = rnorm(NROW(exdata_Upoly)), samples = 3, data = exdata_Upoly, xcoord = xcoord, ycoord = ycoord), "matrix"))
})
