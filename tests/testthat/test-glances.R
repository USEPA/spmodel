load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("glances works geostatistical", {
  spmod1 <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  spmod2 <- splm(y ~ x, exdata, "matern", xcoord, ycoord)
  expect_s3_class(glances(spmod1, spmod2), "tbl")
  expect_equal(NROW(glances(spmod1, spmod2)), 2)
  expect_equal(NCOL(glances(spmod1, spmod2)), 10)
  expect_equal(rbind(glance(spmod1), glance(spmod2)), glances(spmod1, spmod2)[, -1], ignore_attr = TRUE)
})

test_that("glances works autoregressive", {
  spmod1 <- spautor(y ~ x, exdata_poly, "car")
  spmod2 <- spautor(y ~ x, exdata_poly, "sar")
  expect_s3_class(glances(spmod1, spmod2), "tbl")
  expect_equal(NROW(glances(spmod1, spmod2)), 2)
  expect_equal(NCOL(glances(spmod1, spmod2)), 10)
  expect_equal(rbind(glance(spmod2), glance(spmod1)), glances(spmod1, spmod2)[, -1], ignore_attr = TRUE)
})
