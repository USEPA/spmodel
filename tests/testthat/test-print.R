load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("printing works for geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_output(print(spmod))
  expect_output(print(summary(spmod)))
  expect_output(print(anova(spmod)))
})

test_that("printing works for auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_output(print(spmod))
  expect_output(print(summary(spmod)))
  expect_output(print(anova(spmod)))
})
