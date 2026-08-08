test_that("decorrelate works", {

  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spcov_type <- "exponential"
  decorr1 <- decorrelate(y ~ x, exdata, spcov_type = spcov_type, xcoord = xcoord, ycoord = "ycoord")
  expect_type(decorr1, "list")
  preds1 <- predict(decorr1, newdata = newexdata)
  expect_vector(preds1)
  expect_error(tidy(decorr1$grid), NA)

  spcov_type <- "spherical"
  decorr1 <- decorrelate(y ~ x, exdata, spcov_type = spcov_type, xcoord = xcoord, ycoord = "ycoord", local = TRUE)
  expect_type(decorr1, "list")
  preds1 <- predict(decorr1, newdata = newexdata)
  expect_vector(preds1)
  expect_error(tidy(decorr1$grid), NA)

})
