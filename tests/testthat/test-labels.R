load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("labels works geostatistical", {
  spmod <- splm(y ~ x + group, exdata, "exponential", xcoord, ycoord)
  expect_equal(labels(spmod), c("x", "group"))
})

test_that("labels works autoregressive", {
  spmod <- spautor(y ~ x + group, exdata_poly, "car")
  expect_equal(labels(spmod), c("x", "group"))
})
