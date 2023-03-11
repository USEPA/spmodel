################################################################################
############################ check_optim_method (test-check_optim_method.R)
################################################################################

test_that("check optim method works", {

  optim_dotlist <- get_optim_dotlist()

  optim_par <- c(a = 1, b = 2)
  optim_dotlist_val <- check_optim_method(optim_par, optim_dotlist)
  expect_equal(optim_dotlist_val$method, optim_dotlist$method)
  expect_equal(optim_dotlist_val$lower, optim_dotlist$lower)
  expect_equal(optim_dotlist_val$upper, optim_dotlist$upper)

  optim_par <- c(a = 1)
  optim_dotlist_val <- check_optim_method(optim_par, optim_dotlist)
  expect_equal(optim_dotlist_val$method, "Brent")
  expect_equal(optim_dotlist_val$lower, -50)
  expect_equal(optim_dotlist_val$upper, 50)
})

################################################################################
############################ hwInv (test-hwInv.R)
################################################################################

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

################################################################################
############################ randcov_initial (test-randcov_initial.R)
################################################################################

test_that("randcov_initial works", {
  randcov_initial_val <- randcov_initial(group = 1)
  expect_equal(length(randcov_initial_val), 2)
  expect_equal(unname(randcov_initial_val$initial), 1)
  expect_equal(unname(randcov_initial_val$is_known), FALSE)

  randcov_initial_val <- randcov_initial(group = 1, known = "group")
  expect_equal(length(randcov_initial_val), 2)
  expect_equal(unname(randcov_initial_val$initial), 1)
  expect_equal(unname(randcov_initial_val$is_known), TRUE)

  randcov_initial_val <- randcov_initial(1, nm = "group", known = "group")
  expect_equal(length(randcov_initial_val), 2)
  expect_equal(unname(randcov_initial_val$initial), 1)
  expect_equal(unname(randcov_initial_val$is_known), TRUE)

  randcov_initial_val <- randcov_initial(group = 1, known = "given")
  expect_equal(length(randcov_initial_val), 2)
  expect_equal(unname(randcov_initial_val$initial), 1)
  expect_equal(unname(randcov_initial_val$is_known), TRUE)

  randcov_initial_val <- randcov_initial(group = 1, subgroup = 2, known = "group")
  expect_equal(length(randcov_initial_val), 2)
  expect_equal(unname(randcov_initial_val$initial), c(1, 2))
  expect_equal(unname(randcov_initial_val$is_known), c(TRUE, FALSE))

  randcov_initial_val <- randcov_initial(c(1, 2), nm = c("group", "subgroup"), known = "group")
  expect_equal(length(randcov_initial_val), 2)
  expect_equal(unname(randcov_initial_val$initial), c(1, 2))
  expect_equal(unname(randcov_initial_val$is_known), c(TRUE, FALSE))

  randcov_initial_val <- randcov_initial(c(group = 1, subgroup = 2), known = "group")
  expect_equal(length(randcov_initial_val), 2)
  expect_equal(unname(randcov_initial_val$initial), c(1, 2))
  expect_equal(unname(randcov_initial_val$is_known), c(TRUE, FALSE))

  randcov_initial_val <- randcov_initial(c(group = 1, subgroup = 2), known = "xyz")
  expect_equal(length(randcov_initial_val), 2)
  expect_equal(unname(randcov_initial_val$initial), c(1, 2))
  expect_equal(unname(randcov_initial_val$is_known), c(FALSE, FALSE))
})

test_that("errors return", {
  expect_error(randcov_initial(group = NA, known = "given"))
  expect_error(randcov_initial(group = NA, subgroup = 2, known = "given"))
  expect_error(randcov_initial(group = NA, subgroup = 2, known = "subgroup"), NA)
  # expect_error(randcov_initial(c(c(1, 2)))) was for names missing
})

test_that("names of randcov_initial work", {
  randcov_initial1 <- randcov_initial(a = 1, b = 2)
  randcov_initial2 <- randcov_initial(c = 1, d = 2, nm = c("a", "b"))
  randcov_initial3 <- randcov_initial(c(a = 1), b = 2)
  randcov_initial4 <- randcov_initial(c(a = 1), c(b = 2))
  randcov_initial5 <- randcov_initial(1, 2, nm = c("a", "b"))
  expect_equal(randcov_initial1, randcov_initial2)
  expect_equal(randcov_initial1, randcov_initial3)
  expect_equal(randcov_initial1, randcov_initial4)
  expect_equal(randcov_initial1, randcov_initial5)
})

################################################################################
############################ randcov_names (test-randcov_names.R)
################################################################################

test_that("get_randcov_names works", {
  expect_equal(get_randcov_names(~group), "1 | group")
  expect_equal(get_randcov_names(~ 1 | group), "1 | group")
  expect_equal(get_randcov_names(~ group + subgroup), c("1 | group", "1 | subgroup"))
  expect_equal(get_randcov_names(~ group:subgroup), "1 | group:subgroup")
  expect_equal(get_randcov_names(~ group / subgroup), c("1 | group", "1 | group:subgroup"))
  expect_equal(get_randcov_names(~ 1 | group / subgroup), c("1 | group", "1 | group:subgroup"))
  expect_equal(get_randcov_names(~ x | group), c("1 | group", "x | group"))
  expect_equal(get_randcov_names(~ x - 1 | group), c("x | group"))
  expect_equal(get_randcov_names(~ group + (x | group)), c("1 | group", "1 | group", "x | group"))
  expect_equal(get_randcov_names(~ group * x | group), c("1 | group", "group | group", "x | group", "group:x | group"))
  expect_equal(get_randcov_names(~ (x | group) + (x | subgroup)), c("1 | group", "x | group", "1 | subgroup", "x | subgroup"))
  # expect_equal(get_randcov_names(~ group + x | group), c("group + x | group")) # force error here (not intended as () MUST be there)
})

################################################################################
############################ randcov_params (test-randcov_params.R)
################################################################################

test_that("randcov_params works", {
  randcov_params_val <- randcov_params(group = 1, subgroup = 2)
  expect_equal(names(randcov_params_val), c("group", "subgroup"))
  expect_equal(unname(randcov_params_val), c(1, 2))

  randcov_params_val <- randcov_params(1, 2, nm = c("group", "subgroup"))
  expect_equal(names(randcov_params_val), c("group", "subgroup"))
  expect_equal(unname(randcov_params_val), c(1, 2))

  randcov_params_val <- randcov_params(c(1, 2), nm = c("group", "subgroup"))
  expect_equal(names(randcov_params_val), c("group", "subgroup"))
  expect_equal(unname(randcov_params_val), c(1, 2))
})

################################################################################
############################ smwInv_rand (test-smwInv_rand.R)
################################################################################

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
  expect_equal(as.vector(SigInv), as.vector(smwInv_rand_val$SigInv))
  # log det
  expect_equal(as.vector(Sigldet), as.vector(smwInv_rand_val$Sigldet))
})

################################################################################
############################ spcov_initial (test-spcov_initial.R)
################################################################################

test_that("spcov_initial works", {
  spcov_initial_val <- spcov_initial("exponential", de = 1, range = 0.4)
  expect_equal(length(spcov_initial_val), 2)
  expect_equal(unname(spcov_initial_val$initial), c(1, 0.4))
  expect_equal(unname(spcov_initial_val$is_known), c(FALSE, FALSE))

  spcov_initial_val <- spcov_initial("exponential", de = 1, range = 0.4, known = "de")
  expect_equal(length(spcov_initial_val), 2)
  expect_equal(unname(spcov_initial_val$initial), c(1, 0.4))
  expect_equal(unname(spcov_initial_val$is_known), c(TRUE, FALSE))

  spcov_initial_val <- spcov_initial("exponential", de = 1, range = 0.4, known = "given")
  expect_equal(length(spcov_initial_val), 2)
  expect_equal(unname(spcov_initial_val$initial), c(1, 0.4))
  expect_equal(unname(spcov_initial_val$is_known), c(TRUE, TRUE))

  spcov_initial_val <- spcov_initial("exponential", de = 1, range = 0.4, known = "xyz")
  expect_equal(length(spcov_initial_val), 2)
  expect_equal(unname(spcov_initial_val$initial), c(1, 0.4))
  expect_equal(unname(spcov_initial_val$is_known), c(FALSE, FALSE))
})

test_that("names of spcov_initial work", {
  spcov_initial1 <- spcov_initial("matern", de = 1, ie = 1, range = 1, extra = 1, rotate = 1, scale = 1)
  spcov_initial2 <- spcov_initial("matern", de = c(a = 1), ie = c(a = 1), range = c(a = 1), extra = c(a = 1), rotate = c(a = 1), scale = c(a = 1))
  expect_equal(spcov_initial1, spcov_initial2)
})

################################################################################
############################ spcov_matrix (test-spcov_matrix.R)
################################################################################

test_that("spcov_matrix works", {
  h <- spdist(xcoord_val = 1:20)
  spcov1 <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  expect_equal(spcov_matrix(spcov1, h), exp(-h) + diag(nrow(h)))
})

################################################################################
############################ spcov_params (test-spcov_params.R)
################################################################################

test_that("spcov_params works", {
  spcov1 <- spcov_params("exponential", de = 1, ie = 2, range = 1, rotate = 2, scale = 0.75)
  expect_s3_class(spcov1, "exponential")
  expect_equal(length(spcov1), 5)
  expect_equal(names(spcov1), c("de", "ie", "range", "rotate", "scale"))
  spcov2 <- spcov_params("matern", de = 1, ie = 2, range = 1, extra = 1 / 2, rotate = 1, scale = 0.5)
  expect_s3_class(spcov2, "matern")
  expect_equal(length(spcov2), 6)
  expect_equal(names(spcov2), c("de", "ie", "range", "extra", "rotate", "scale"))
})

test_that("spcov_params works", {
  expect_error(spcov_params("exponential", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("spherical", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("gaussian", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("triangular", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("circular", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("cubic", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("pentaspherical", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("cosine", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("wave", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("jbessel", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("gravity", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("rquad", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("magnetic", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("matern", de = 1, ie = 1, range = 1, extra = 1), NA)
  expect_error(spcov_params("cauchy", de = 1, ie = 1, range = 1, extra = 1), NA)
  expect_error(spcov_params("pexponential", de = 1, ie = 1, range = 1, extra = 1), NA)
  expect_error(spcov_params("car", de = 1, ie = 1, range = 1, extra = 1), NA)
  expect_error(spcov_params("sar", de = 1, ie = 1, range = 1, extra = 1), NA)
  expect_error(spcov_params("none", ie = 1), NA)
})

test_that("spcov_params errors", {
  expect_error(spcov_params("exponential", de = -1, ie = 1, range = 1))
  expect_error(spcov_params("exponential", de = 1, ie = -1, range = 1))
  expect_error(spcov_params("exponential", de = 1, ie = -1, range = -1))
  expect_error(spcov_params("exponential", de = 1, ie = 1, range = 1, rotate = 10))
  expect_error(spcov_params("exponential", de = 1, ie = 1, range = 1, rotate = 1, scale = 2))
  expect_error(spcov_params("exponential", de = 1, ie = 1, range = 1, rotate = -10))
  expect_error(spcov_params("exponential", de = 1, ie = 1, range = 1, rotate = 1, scale = -2))
  expect_error(spcov_params("matern", de = 1, ie = -1, range = 1, extra = -1))
  expect_error(spcov_params("matern", de = 1, ie = -1, range = 1, extra = 1 / 6))
  expect_error(spcov_params("matern", de = 1, ie = -1, range = 1, extra = 6))
  expect_error(spcov_params("cauchy", de = 1, ie = 1, range = 1, extra = -5))
  expect_error(spcov_params("cauchy", de = 1, ie = 1, range = 1, extra = 0))
  expect_error(spcov_params("pexponential", de = 1, ie = 1, range = 1, extra = 3))
  expect_error(spcov_params("pexponential", de = 1, ie = 1, range = 1, extra = 0))
  expect_error(spcov_params("pexponential", de = 1, ie = 1, range = 1, extra = -3))
  expect_error(spcov_params("car", de = 1, range = 1, ie = 1))

  # spcov type problems
  expect_error(spcov_params(de = 1, ie = 1, range = 1))
  expect_error(spcov_params("xyz", de = 1, ie = 1, range = 1))
})

test_that("defaults are applied", {
  expect_error(spcov_params("exponential", de = 1, ie = 1, range = 1), NA)
  expect_error(spcov_params("car", de = 1, range = 1, extra = 1), NA)
})

################################################################################
############################ spcov_vector (test-spcov_vector.R)
################################################################################

test_that("spcov_vector works", {
  h <- 1:20
  spcov1 <- spcov_params("exponential", de = 1, ie = 1, range = 1)
  expect_equal(spcov_vector.exponential(spcov1, h), exp(-h))
})
