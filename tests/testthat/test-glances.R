load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_poly.rda", package = "spmodel"))


test_local <- TRUE # FALSE for CRAN

##### CRAN test
test_that("glances works geostatistical", {
  spmod1 <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  spmod2 <- splm(y ~ x, exdata, "none", xcoord, ycoord)
  expect_s3_class(glances(spmod1, spmod2), "tbl")
  expect_equal(NROW(glances(spmod1, spmod2)), 2)
  expect_equal(NCOL(glances(spmod1, spmod2)), 10)
  expect_equal(rbind(glance(spmod1), glance(spmod2)), glances(spmod1, spmod2)[, -1], ignore_attr = TRUE)
})

#### local tests
if (test_local) {
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

  test_that("glances works iteratively splm spcov_type", {
    spmod <- splm(y ~ x, exdata, c("exponential", "matern"), xcoord, ycoord)
    expect_s3_class(glances(spmod), "tbl")
    expect_equal(NROW(glances(spmod)), 2)
    expect_equal(NCOL(glances(spmod)), 10)
    expect_equal(rbind(glance(spmod$exponential), glance(spmod$matern)), glances(spmod)[, -1], ignore_attr = TRUE)
  })

  test_that("glances works iteratively splm spcov_initial", {
    spcov_init <- lapply(c("exponential", "matern"), function(x) spcov_initial(x, de = 1))
    spmod <- splm(y ~ x, exdata, spcov_initial = spcov_init, xcoord = xcoord, ycoord = ycoord)
    expect_s3_class(glances(spmod), "tbl")
    expect_equal(NROW(glances(spmod)), 2)
    expect_equal(NCOL(glances(spmod)), 10)
    expect_equal(rbind(glance(spmod$spcov_initial_1), glance(spmod$spcov_initial_2)), glances(spmod)[, -1], ignore_attr = TRUE)
  })

  test_that("glances works iteratively spautor spcov_type", {
    spmod <- spautor(y ~ x, exdata_poly, c("car", "sar"))
    expect_s3_class(glances(spmod), "tbl")
    expect_equal(NROW(glances(spmod)), 2)
    expect_equal(NCOL(glances(spmod)), 10)
    expect_equal(rbind(glance(spmod$car), glance(spmod$sar)), glances(spmod, sort_by = "order")[, -1], ignore_attr = TRUE)
    # sort by order here because sort by is off in natural implementation by AICc
  })

  test_that("glances works iteratively spautor spcov_initial", {
    spcov_init <- lapply(c("car", "sar"), function(x) spcov_initial(x, de = 1))
    spmod <- spautor(y ~ x, exdata_poly, spcov_initial = spcov_init)
    expect_s3_class(glances(spmod), "tbl")
    expect_equal(NROW(glances(spmod)), 2)
    expect_equal(NCOL(glances(spmod)), 10)
    expect_equal(rbind(glance(spmod$spcov_initial_1), glance(spmod$spcov_initial_2)), glances(spmod, sort_by = "order")[, -1], ignore_attr = TRUE)
  })

}
