load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("glance works geostatistical", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_s3_class(glance(spmod), "tbl")
  expect_equal(NROW(glance(spmod)), 1)
  expect_equal(NCOL(glance(spmod)), 9)
  expect_false(any(is.na(glance(spmod))))

  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, estmethod = "sv-wls")
  expect_s3_class(glance(spmod), "tbl")
  expect_equal(NROW(glance(spmod)), 1)
  expect_equal(NCOL(glance(spmod)), 9)
  expect_true(any(is.na(glance(spmod))))
})

test_that("glance works autoregressive", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_s3_class(glance(spmod), "tbl")
  expect_equal(NCOL(glance(spmod)), 9)
  expect_false(any(is.na(glance(spmod))))
})
