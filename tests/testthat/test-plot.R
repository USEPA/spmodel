load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))

test_that("plot works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)

  # plot 1
  expect_error(plot(spmod, which = 1), NA)

  # plot 2
  expect_error(plot(spmod, which = 2), NA)

  # plot 3
  expect_error(plot(spmod, which = 3), NA)

  # plot 4
  expect_error(plot(spmod, which = 4), NA)

  # plot 5
  expect_error(plot(spmod, which = 5), NA)

  # plot 6
  expect_error(plot(spmod, which = 6), NA)

  # plot 7
  expect_error(plot(spmod, which = 7), NA)

  # plot 8
  expect_error(plot(spmod, which = 8), NA)

  # plot 9 (return error)
  expect_error(plot(spmod, which = 9))

  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord, anisotropy = TRUE)

  # plot 8
  expect_error(plot(spmod, which = 8), NA)
})

test_that("plot works auto", {
  spmod <- spautor(y ~ x, exdata_poly, "car")

  # plot 1
  expect_error(plot(spmod, which = 1), NA)

  # plot 2
  expect_error(plot(spmod, which = 2), NA)

  # plot 3
  expect_error(plot(spmod, which = 3), NA)

  # plot 4
  expect_error(plot(spmod, which = 4), NA)

  # plot 5
  expect_error(plot(spmod, which = 5), NA)

  # plot 6
  expect_error(plot(spmod, which = 6), NA)

  # plot 7
  expect_error(plot(spmod, which = 7))

  # plot 8
  expect_error(plot(spmod, which = 8))
})
