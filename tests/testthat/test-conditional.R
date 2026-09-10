test_that("conditional works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  cond1 <- conditional(spmod, newdata = newexdata, samples = 50)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))

  cond_all <- conditional(spmod, newdata = newexdata, output = "all", samples = 50)
  expect_type(cond_all, "list")
  expect_named(cond_all, c("newdata", "beta", "object"))
  expect_equal(dim(cond_all$newdata), c(NROW(newexdata), 50))
  expect_equal(dim(cond_all$beta), c(length(coef(spmod)), 50))
  expect_equal(dim(cond_all$object), c(spmod$n, 50))

  # falls back to object$newdata (the missing rows) when newdata is omitted
  exdata_miss <- exdata
  exdata_miss$y[1:5] <- NA
  spmod_miss <- splm(y ~ x, exdata_miss, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  cond_miss <- conditional(spmod_miss, samples = 50)
  expect_equal(nrow(cond_miss), 5)
})

test_that("conditional works for spglm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  exdata_pois <- exdata
  exdata_pois$count <- round(abs(exdata_pois$y) * 3)
  spmod <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  cond_link <- conditional(spmod, newdata = newexdata, samples = 50)
  expect_equal(dim(cond_link), c(NROW(newexdata), 50))

  cond_response <- conditional(spmod, newdata = newexdata, type = "response", samples = 50)
  expect_true(all(cond_response >= 0))

  cond_new <- conditional(spmod, newdata = newexdata, type = "new", samples = 50)
  expect_true(all(cond_new >= 0))
  expect_equal(cond_new, round(cond_new)) # poisson draws are counts
})

test_that("conditional() simulate_covparams = TRUE works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_false(is.null(spmod$vcov$cov)) # ddf = "satterthwaite" default for n <= 500

  cond1 <- conditional(spmod, newdata = newexdata, samples = 50, simulate_covparams = TRUE)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))
})

test_that("conditional() output = 'cov'/'spcov'/'randcov' works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  spcov_names <- names(coef(spmod, type = "spcov"))

  cov_out <- conditional(spmod, newdata = newexdata, samples = 30, simulate_covparams = TRUE, output = "cov")
  expect_equal(dim(cov_out), c(length(spcov_names), 30))
  expect_equal(rownames(cov_out), spcov_names)
  expect_true(all(is.finite(cov_out)))
})
