load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("vcov works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_true(isSymmetric(vcov(spmod)))
  expect_true(all(diag(vcov(spmod)) >= 0))
})

test_that("vcov works auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_true(isSymmetric(vcov(spmod)))
  expect_true(all(diag(vcov(spmod)) >= 0))
})
