load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("hatvalues works geostatistical", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_vector(hatvalues(spmod))
  expect_true(min(hatvalues(spmod)) >= 0)
  expect_true(max(hatvalues(spmod)) <= 1)
  expect_equal(sum(hatvalues(spmod)), spmod$p)
})

test_that("hatvalues works autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_vector(hatvalues(spmod))
  expect_true(min(hatvalues(spmod)) >= 0)
  expect_true(max(hatvalues(spmod)) <= 1)
  expect_equal(sum(hatvalues(spmod)), spmod$p)
})
