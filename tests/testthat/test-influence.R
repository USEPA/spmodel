load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("influence works geostatistical", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_error(influence(spmod), NA)
  expect_equal(NROW(influence(spmod)), 100)
  expect_equal(NCOL(influence(spmod)), 4)
  expect_equal(colnames(influence(spmod)), c(".resid", ".hat", ".cooksd", ".std.resid"))
  expect_s3_class(influence(spmod), "tbl")
})

test_that("influence works autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_error(influence(spmod), NA)
  expect_equal(NROW(influence(spmod)), 49)
  expect_equal(NCOL(influence(spmod)), 4)
  expect_equal(colnames(influence(spmod)), c(".resid", ".hat", ".cooksd", ".std.resid"))
  expect_s3_class(influence(spmod), "tbl")
})
