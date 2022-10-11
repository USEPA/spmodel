load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

test_local <- FALSE # FALSE for CRAN


if (!requireNamespace("ranger", quietly = TRUE)) {
  expect_equal(2, 2) # dummy test
} else {

  #### CRAN check

  test_that("the model runs", {
    spcov_type <- "car"
    expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type), NA)
  })

  test_that("prediction works", {
    spcov_type <- "car"
    sprfmod <- spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type)
    expect_vector(predict(sprfmod))
  })


}
