test_that("hwInv works", {
  h <- as.matrix(dist(1:10))
  Sig <- exp(-h / 10)
  SigInv <- solve(Sig)
  Sigldet <- as.numeric(determinant(Sig)$modulus)

  observed_index <- 3:10
  Sig_small <- Sig[observed_index, observed_index]
  SigInv_small <- solve(Sig_small)
  Sigldet_small <- as.numeric(determinant(Sig_small)$modulus)

  hwInv_val <- hwInv(SigInv, Sigldet, observed_index)

  expect_equal(unname(round(as.matrix(Sig_small %*% SigInv_small), 12)), diag(8))
  expect_equal(unname(round(as.matrix(Sig_small %*% hwInv_val$SigInv), 12)), diag(8)) # some floating point makes Matrix version of inverse not
  # exactly equal to base version
  expect_equal(Sigldet_small, hwInv_val$Sigldet)
})
