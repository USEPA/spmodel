load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

test_local <- FALSE # FALSE for CRAN

#### CRAN checks
test_that("covmatrix works geostatistical", {
  spcov_type <- "exponential"
  spmod <- splm(y ~ x, exdata_M, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord)
  expect_equal(c(99, 99), dim(covmatrix(spmod)))
  expect_equal(c(1, 99), dim(covmatrix(spmod, newdata = spmod$newdata)))
  expect_equal(c(10, 99), dim(covmatrix(spmod, newdata = newexdata)))
})

test_that("covmatrix works autoregressive", {
  spcov_type <- "car"
  spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type)
  expect_equal(c(48, 48), dim(covmatrix(spmod)))
  expect_equal(c(1, 48), dim(covmatrix(spmod, newdata = spmod$newdata)))
})


if (test_local) {

  #### CRAN checks
  test_that("covmatrix works geostatistical", {
    spcov_type <- "exponential"
    spmod <- splm(y ~ x, exdata_M, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, random = ~ group)
    expect_equal(c(99, 99), dim(covmatrix(spmod)))
    expect_equal(c(1, 99), dim(covmatrix(spmod, newdata = spmod$newdata)))
    expect_equal(c(10, 99), dim(covmatrix(spmod, newdata = newexdata)))
  })

  test_that("covmatrix works autoregressive", {
    spcov_type <- "car"
    spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type, random = ~ group)
    expect_equal(c(48, 48), dim(covmatrix(spmod)))
    expect_equal(c(1, 48), dim(covmatrix(spmod, newdata = spmod$newdata)))
  })

  test_that("covmatrix works geostatistical", {
    spcov_type <- "exponential"
    spmod <- splm(y ~ x, exdata_M, spcov_type = spcov_type, xcoord = xcoord, ycoord = ycoord, partition_factor = ~ group)
    expect_equal(c(99, 99), dim(covmatrix(spmod)))
    expect_equal(c(1, 99), dim(covmatrix(spmod, newdata = spmod$newdata)))
    expect_equal(c(10, 99), dim(covmatrix(spmod, newdata = newexdata)))
  })

  test_that("covmatrix works autoregressive", {
    spcov_type <- "car"
    spmod <- spautor(y ~ x, exdata_Mpoly, spcov_type = spcov_type, partition_factor = ~ group)
    expect_equal(c(48, 48), dim(covmatrix(spmod)))
    expect_equal(c(1, 48), dim(covmatrix(spmod, newdata = spmod$newdata)))
  })

}
