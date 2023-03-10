load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("residuals works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_vector(residuals(spmod))
  expect_equal(residuals(spmod), exdata$y - fitted(spmod))
  expect_vector(residuals(spmod, type = "response"))
  expect_vector(residuals(spmod, type = "pearson"))
  expect_vector(residuals(spmod, type = "standardized"))
  expect_equal(residuals(spmod, type = "standardized"), rstandard(spmod))
})

test_that("residuals works auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_vector(residuals(spmod))
  expect_equal(residuals(spmod), exdata_poly$y - fitted(spmod))
  expect_vector(residuals(spmod, type = "response"))
  expect_vector(residuals(spmod, type = "pearson"))
  expect_vector(residuals(spmod, type = "standardized"))
  expect_equal(residuals(spmod, type = "standardized"), rstandard(spmod))
})
