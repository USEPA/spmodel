skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))
skip_if_not_installed("GPvecchia")
skip_if_not_installed("spsurvey")

coincident_fixture <- function(glm = FALSE, observed = TRUE, prediction = TRUE, random = FALSE) {
  set.seed(1)
  d <- data.frame(cx = if (observed) rep(0:5, each = 2) else seq(0, 5, length.out = 12),
    cy = 0, x = rep(c(-1, 1), 6), g = factor(rep(1:3, 4)),
    y = c(2, 1, 3, 4, 2, 5, 1, 2, 4, 3, 1, 2))
  new <- d[c(1, 2, 7), ]
  new$cx <- if (prediction) c(.1, .1, 3.1) else c(.1, 1.1, 3.1)
  args <- list(formula = y ~ x, data = d, xcoord = "cx", ycoord = "cy",
    spcov_initial = spcov_initial("exponential", de = .4, ie = .2, range = .15, known = "given"))
  if (glm) args$family <- "poisson"
  if (random) {
    args$random <- ~ (1 | g)
    args$randcov_initial <- randcov_initial(g = .05, known = "given")
  }
  list(fit = do.call(if (glm) spglm else splm, args), new = new, data = d)
}

test_that("coincident simulation defaults are finite and retain covariance", {
  set.seed(1)
  d <- data.frame(cx = rep(0:5, each = 2), cy = 0)
  params <- spcov_params("exponential", de = 1, ie = .2, range = .2)
  internal <- setNames(d, c("...xcoord...", "...ycoord..."))
  settings <- get_local_list_simulation_vecchia(list(), nrow(d), internal)
  expect_identical(settings$ordering, "grts")
  expect_equal(sort(settings$order), seq_len(nrow(d)))
  set.seed(1)
  draws <- sprnorm(params, data = d, xcoord = cx, ycoord = cy, samples = 6000, local = TRUE)
  expect_true(all(is.finite(draws)))
  expect_lt(max(abs(apply(draws, 1, var) / 1.2 - 1)), .07)
  expect_lt(abs(cor(draws[1, ], draws[2, ]) - 1 / 1.2), .025)
  set.seed(1)
  explicit <- sprnorm(params, data = d, xcoord = cx, ycoord = cy, samples = 6000,
    local = list(ordering = "grts"))
  expect_identical(draws, explicit)
  expect_true(all(is.finite(sprpois(params, data = d, xcoord = cx, ycoord = cy,
    samples = 3, local = TRUE))))
  expect_error(sprnorm(params, data = d, xcoord = cx, ycoord = cy,
    local = list(ordering = "maxmin")), "coincident")
  expect_error(sprpois(params, data = d, xcoord = cx, ycoord = cy,
    local = list(ordering = "maxmin")), "coincident")
  expect_equal(sort(conditional_vecchia_order(matrix(0, 3, 2))), 1:3)
})

test_that("decorrelation and block prediction reject explicit coincident maxmin", {
  f <- coincident_fixture()
  params <- coef(f$fit, type = "spcov")
  expect_error(decorrelate_data(y ~ x, f$data, params, xcoord = cx, ycoord = cy,
    ordering = "maxmin"), "coincident")
  expect_error(decorrelate(y ~ x, f$data, spcov_params = params, xcoord = cx, ycoord = cy,
    ordering = "maxmin"), "coincident")
  expect_error(predict(f$fit, f$new, block = TRUE,
    local = list(method_new = "subset", size_new = 2, ordering = "maxmin")), "coincident")
  set.seed(1)
  default <- decorrelate_data(y ~ x, f$data, params, xcoord = cx, ycoord = cy)
  set.seed(1)
  explicit <- decorrelate_data(y ~ x, f$data, params, xcoord = cx, ycoord = cy, ordering = "grts")
  expect_identical(default, explicit)
})

test_that("conditional defaults resolve each set independently and preserve empty output", {
  for (observed in c(FALSE, TRUE)) for (prediction in c(FALSE, TRUE)) {
    f <- coincident_fixture(glm = TRUE, observed = observed, prediction = prediction, random = TRUE)
    set.seed(1)
    local <- get_local_list_conditional(list(size = 5), f$fit, f$new)
    expect_identical(local$ordering, if (prediction) "grts" else "maxmin")
    expect_identical(local$ordering_observed, if (observed) "grts" else "maxmin")
    if (observed || prediction) expect_error(conditional(f$fit, f$new, samples = 2,
      local = list(size = 5, ordering = "maxmin")), "coincident")
    set.seed(1)
    no_new <- conditional(f$fit, samples = 2, output = c("latent", "random"), local = list(size = 5))
    expect_identical(rownames(no_new$latent), rownames(model.matrix(f$fit)))
    expect_equal(dim(no_new$latent), c(12, 2))
    expect_true(all(is.finite(no_new$latent)))
    expect_true(all(is.finite(no_new$random)))
    one_new <- conditional(f$fit, f$new[1, , drop = FALSE], samples = 2, local = list(size = 5))
    expect_equal(dim(one_new), c(1, 2))
    expect_true(all(is.finite(one_new)))
  }
  # LM observed coordinates are not ordered; only its distinct prediction set is.
  f <- coincident_fixture(prediction = FALSE)
  expect_true(all(is.finite(conditional(f$fit, f$new, samples = 2,
    local = list(size = 5, ordering = "maxmin")))))
})

test_that("coincident conditional moments agree with exact", {
  for (glm in c(FALSE, TRUE)) for (scenario in c("observed", "prediction", "both", "random")) {
    f <- coincident_fixture(glm = glm, observed = scenario != "prediction",
      prediction = scenario != "observed", random = scenario == "random")
    output <- c("newdata", "beta", if (glm) "latent", if (scenario == "random") "random")
    set.seed(1)
    exact <- conditional(f$fit, f$new, output = output, samples = 6000, local = FALSE)
    set.seed(1)
    local <- conditional(f$fit, f$new, output = output, samples = 6000, local = list(size = 5))
    expect_identical(rownames(local$newdata), rownames(exact$newdata))
    if (glm) expect_identical(rownames(local$latent), rownames(model.matrix(f$fit)))
    a <- do.call(rbind, exact); b <- do.call(rbind, local)
    expect_true(all(is.finite(b)))
    s <- apply(a, 1, sd)
    # Includes Monte Carlo differences and finite-neighborhood approximation.
    expect_lt(max(abs(rowMeans(b) - rowMeans(a)) / s), .10)
    expect_lt(max(abs(apply(b, 1, sd) / s - 1)), .12)
    expect_lt(max(abs(cov(t(b)) - cov(t(a))) / outer(s, s)), .15)
    set.seed(1)
    all <- conditional(f$fit, f$new, output = output, samples = 3,
      local = list(method = "all", ordering = "maxmin"))
    set.seed(1)
    exact_all <- conditional(f$fit, f$new, output = output, samples = 3, local = FALSE)
    expect_identical(all, exact_all)
  }
})
