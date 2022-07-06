load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("free generics work geostatistical", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_true(is.language(terms(spmod))) # works because there is an object spmod$terms
  expect_true(is.call(getCall(spmod))) # works because there is an object spmod$call
  expect_equal(update(formula(spmod), y ~ 1), y ~ 1) # works because formula(object) works (object has class spmod)
  spmod <- update(spmod, y ~ 1, spcov_type = "spherical")
  expect_s3_class(spmod, "spmod")
  expect_equal(spmod$fn, "splm")
  expect_equal(formula(spmod), y ~ 1)
  expect_s3_class(coefficients(spmod, type = "spcov"), "spherical")
  spmod <- update(spmod, . ~ . + offset(x))
  expect_vector(model.offset(model.frame(spmod))) # works because model.frame(object) works (object has class spmod)
  expect_vector(model.response(model.frame(spmod))) # works because model.frame(object) works (object has class spmod)
})

test_that("free generics work autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_true(is.language(terms(spmod))) # works because there is an object spmod$terms
  expect_true(is.call(getCall(spmod))) # works because there is an object spmod$call
  expect_equal(update(formula(spmod), y ~ 1), y ~ 1) # works because formula(object) works (object has class spmod)
  spmod <- update(spmod, y ~ 1, spcov_type = "sar")
  expect_s3_class(spmod, "spmod")
  expect_equal(spmod$fn, "spautor")
  expect_equal(formula(spmod), y ~ 1)
  expect_s3_class(coefficients(spmod, type = "spcov"), "sar")
  spmod <- update(spmod, . ~ . + offset(x))
  expect_vector(model.offset(model.frame(spmod))) # works because model.frame(object) works (object has class spmod)
  expect_vector(model.response(model.frame(spmod))) # works because model.frame(object) works (object has class spmod)
})
