load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("model.frame works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_equal(NCOL(model.frame(spmod)), 2)
  expect_equal(NROW(model.frame(spmod)), 100)
  expect_equal(model.frame(spmod), model.frame(y ~ x, exdata, drop.unused.levels = TRUE, na.action = na.omit))
})

test_that("model.frame works auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_equal(NCOL(model.frame(spmod)), 2)
  expect_equal(NROW(model.frame(spmod)), 49)
  expect_equal(model.frame(spmod), model.frame(y ~ x, exdata_poly, drop.unused.levels = TRUE, na.action = na.omit))
})
