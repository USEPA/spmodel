load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("pseudoR2 works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_vector(pseudoR2(spmod))
  expect_vector(pseudoR2(spmod, adjust = TRUE))
})

test_that("pseudoR2 works auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")
  expect_vector(pseudoR2(spmod))
  expect_vector(pseudoR2(spmod, adjust = TRUE))
})

test_that("pseudoR2 matches classical r-squared", {
  lmod <- lm(y ~ x, exdata)
  spmod <- splm(y ~ x, exdata, "none")
  expect_equal(summary(lmod)$r.squared, pseudoR2(spmod))
  expect_equal(summary(lmod)$adj.r.squared, pseudoR2(spmod, adjust = TRUE))
})
