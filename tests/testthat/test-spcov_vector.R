test_that("spcov_vector works", {
  h <- 1:20
  spcov1 <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  expect_equal(spcov_vector.exponential(spcov1, h), exp(-h))
})
