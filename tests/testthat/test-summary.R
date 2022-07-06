load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("summary works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_error(summary(spmod), NA)
})

test_that("summary works auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_error(summary(spmod), NA)
})
