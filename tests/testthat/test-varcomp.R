load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))
load(system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

test_local <- FALSE

#### CRAN test
test_that("varcomp works geostatistical", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, estmethod = "reml")
  varcomp_val <- varcomp(spmod)
  expect_equal(NROW(varcomp_val), 3)
})


#### local tests
if (test_local) {

  test_that("varcomp works geostatistical random", {
    spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord,
                  ycoord = ycoord, estmethod = "reml", random = ~ group)
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val), 4)
  })

  test_that("varcomp works autoregressive", {
    spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val), 3)
  })

  test_that("varcomp works autoregressive random", {
    spmod <- spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml", random = ~ group)
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val), 4)
  })

  test_that("varcomp works unconnected autoregressive", {
    spmod <- spautor(y ~ x, exdata_Upoly, spcov_type = "car", estmethod = "reml")
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val$connected), 3)
    expect_equal(NROW(varcomp_val$unconnected), 2)
    expect_equal(length(varcomp_val$ratio), 1)
  })

  test_that("varcomp works unconnected autoregressive random", {
    spmod <- spautor(y ~ x, exdata_Upoly, spcov_type = "car", estmethod = "reml", random = ~ group)
    varcomp_val <- varcomp(spmod)
    expect_equal(NROW(varcomp_val$connected), 4)
    expect_equal(NROW(varcomp_val$unconnected), 3)
    expect_equal(length(varcomp_val$ratio), 1)
  })

}
