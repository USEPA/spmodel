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
