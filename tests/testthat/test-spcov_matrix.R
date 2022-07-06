test_that("spcov_matrix works", {
  h <- spdist(xcoord_val = 1:20)
  spcov1 <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  expect_equal(spcov_matrix(spcov1, h), exp(-h) + diag(nrow(h)))
})
