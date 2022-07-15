load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("sp object error tests geostatistical", {
  # first we make "fake" sp objects that can be loaded without also loading sp
  exdata_sp <- exdata
  attr(class(exdata_sp), "package") <- "sp"
  expect_error(splm(y ~ x, exdata_sp, "exponential", xcoord = xcoord, ycoord = ycoord))
  expect_error(esv(y ~ x, exdata_sp, xcoord = xcoord, ycoord = ycoord))
  spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  expect_error(sprnorm(spcov_params_val, data = exdata_sp, xcoord = xcoord, ycoord = ycoord))

  spmod <- splm(y ~ x, exdata, "exponential", xcoord = xcoord, ycoord = ycoord)
  newexdata_sp <- newexdata
  attr(class(newexdata_sp), "package") <- "sp"
  expect_error(predict(spmod, newexdata_sp))
  expect_error(augment(spmod, newdata = newexdata_sp))
})

test_that("sp object error tests autoregressive", {
  # first we make "fake" sp objects that can be loaded without also loading sp
  exdata_poly_sp <- exdata_poly
  attr(class(exdata_poly_sp), "package") <- "sp"
  expect_error(spautor(y ~ x, exdata_poly_sp, "car"))
  spcov_params_val <- spcov_params("car", de = 1, ie = 1, range = 0.5, extra = 0)
  expect_error(sprnorm(spcov_params_val, data = exdata_poly_sp))
})
