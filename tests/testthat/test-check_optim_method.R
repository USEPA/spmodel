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
