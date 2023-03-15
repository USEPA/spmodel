test_that("blank test", {
  expect_null(NULL)
})

test_local <- FALSE # FALSE for CRAN

if (test_local) {

  set.seed(1)

  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

  if (!requireNamespace("ranger", quietly = TRUE)) {
    expect_equal(2, 2) # dummy test
  } else {

    #### CRAN check

    test_that("the model runs", {
      spcov_type <- "car"
      expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type), NA)
      expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "ml"), NA)
    })

    if (test_local) {
      test_that("the model runs", {
        spcov_type <- "car"
        num.tree <- 499
        spcov_initial_val <- spcov_initial(spcov_type = spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
        expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val), NA)
        expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
      })

      test_that("the model list runs", {
        spcov_type <- c("car", "sar")
        expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type), NA)
        expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_type = spcov_type, estmethod = "ml"), NA)
      })

      test_that("the model list runs", {
        spcov_type <- c("car", "sar")
        num.tree <- 499
        spcov_initial_val <- lapply(spcov_type, function(x) spcov_initial(spcov_type = x, de = 1, ie = 1, range = 0.5, known = "de"))
        expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val), NA)
        expect_error(spautorRF(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
      })
    }

  }

}
