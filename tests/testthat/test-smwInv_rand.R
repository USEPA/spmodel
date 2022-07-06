test_that("smwInv_rand works", {
  x <- 1:32 # x must be even
  h <- as.matrix(dist(x))
  A <- exp(-h)
  data <- data.frame(B1 = factor(rep(c(1, 0), each = length(x) / 2)), B2 = factor(rep(c(1, 0), times = length(x) / 2)))
  randcov_names_val <- get_randcov_names(~ B1 + B2)
  randcov_params_val <- randcov_params(1.5, 2)
  names(randcov_params_val) <- randcov_names_val
  randcov_Zs <- get_randcov_Zs(data, randcov_names_val)

  Ainv <- solve(A)
  Aldet <- as.numeric(determinant(A)$modulus)

  Sig <- A + randcov_params_val[1] * randcov_Zs$`1 | B1`$ZZt + randcov_params_val[2] * randcov_Zs$`1 | B2`$ZZt
  SigInv <- solve(Sig)
  Sigldet <- as.numeric(determinant(Sig)$modulus)

  smwInv_rand_val <- smwInv_rand(Ainv, Aldet, randcov_params = randcov_params_val, randcov_Zs = randcov_Zs)

  # inv
  expect_equal(SigInv, unname(smwInv_rand_val$SigInv))
  # log det
  expect_equal(Sigldet, unname(smwInv_rand_val$Sigldet))
})
