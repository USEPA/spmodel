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
