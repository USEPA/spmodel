# Snapshot tests capturing the current print()/print(summary()) output for
# splm/spautor/spglm/spgautor

test_that("splm print output is unchanged (point-referenced, default)", {
  spmod <- splm(z ~ water + tarp,
    data = caribou,
    spcov_type = "exponential", xcoord = x, ycoord = y
  )
  expect_snapshot(print(spmod))
  expect_snapshot(print(summary(spmod)))
})

test_that("splm print output is unchanged (anisotropy shown)", {
  spmod <- splm(z ~ water + tarp,
    data = caribou,
    spcov_type = "exponential", xcoord = x, ycoord = y,
    anisotropy = TRUE
  )
  expect_snapshot(print(spmod))
})

test_that("splm print output is unchanged (none/ie collapse)", {
  spmod <- splm(z ~ water + tarp,
    data = caribou,
    spcov_type = "none"
  )
  expect_snapshot(print(spmod))
})

test_that("spautor print output is unchanged (de/range only)", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")
  expect_snapshot(print(spmod))
  expect_snapshot(print(summary(spmod)))
})

test_that("spautor print output is unchanged (de/ie/range/extra all present)", {
  load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))
  spmod <- spautor(y ~ x, exdata_Upoly, spcov_type = "car", estmethod = "reml")
  expect_snapshot(print(spmod))
  expect_snapshot(print(summary(spmod)))
})

test_that("spglm print output is unchanged", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)
  spmod <- spglm(y_pois ~ x,
    family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord,
    spcov_type = "exponential"
  )
  expect_snapshot(print(spmod))
  expect_snapshot(print(summary(spmod)))
})

test_that("spgautor print output is unchanged (de/ie/range/extra all present)", {
  load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))
  spmod <- spgautor(abs(y) ~ x, "Gamma", exdata_Upoly, spcov_type = "car", estmethod = "reml")
  expect_snapshot(print(spmod))
  expect_snapshot(print(summary(spmod)))
})
