skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

test_that("one-dimensional simulation warns only for supplied y-coordinates", {
  d <- data.frame(cx = seq(0, 1, length.out = 12), cy = 0)
  settings <- list(FALSE, list(approximation = "vecchia", size = 4, ordering = "none"),
    list(approximation = "low-rank", size_base = 6, size_new = 3,
      reorder_base = "random", kmeans_new = FALSE))
  for (type in c("triangular", "cosine", "exponential")) {
    params <- spcov_params(type, de = .4, ie = .2, range = .3)
    for (local in settings) for (sim in list(sprnorm, sprpois, sprbinom, sprnbinom,
      sprbeta, sprgamma, sprinvgauss)) {
      set.seed(2)
      expect_warning(omitted <- sim(params, data = d, xcoord = cx,
        samples = 3, local = local), NA)
      set.seed(2)
      expect_warning(null <- sim(params, data = d, xcoord = cx, ycoord = NULL,
        samples = 3, local = local), NA)
      expect_identical(null, omitted)
      set.seed(2)
      warning <- if (type == "exponential") NA else "Ignoring y-coordinate"
      expect_warning(supplied <- sim(params, data = d, xcoord = cx, ycoord = cy,
        samples = 3, local = local), warning)
      expect_identical(supplied, omitted)
      expect_true(all(is.finite(omitted)))
    }
  }
})

test_that("low-rank sf simulation preserves one-dimensional covariance", {
  set.seed(3)
  d <- data.frame(cx = runif(40), cy = runif(40))
  dsf <- sf::st_as_sf(d, coords = c("cx", "cy"))
  for (type in c("triangular", "cosine")) for (seed in 1:3) {
    params <- spcov_params(type, de = 1, ie = .3, range = .3)
    for (sim in list(sprnorm, sprpois)) for (kmeans in c(FALSE, TRUE)) {
      local <- list(approximation = "low-rank", size_base = 10, size_new = 5,
        reorder_base = "random", kmeans_new = kmeans)
      set.seed(seed)
      expect_warning(actual <- sim(params, data = dsf, samples = 3, local = local),
        "Ignoring y-coordinate")
      expect_true(all(is.finite(actual)))
      if (!kmeans) {
        set.seed(seed)
        expected <- sim(params, data = d, xcoord = cx, samples = 3, local = local)
        expect_equal(actual, expected, tolerance = 1e-12)
      }
    }
  }
})

test_that("sprnorm() local$approximation = 'vecchia' works", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  spcov_params_val <- spcov_params("exponential", de = 1, ie = 0.05, range = 1.5)

  sim1 <- sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "vecchia", size = 10))
  expect_true(inherits(sim1, "matrix"))
  expect_equal(dim(sim1), c(NROW(exdata), 20))
  expect_true(all(is.finite(sim1)))

  # method = "all" (no truncation) should closely match the exact (local =
  # FALSE) unconditional distribution -- see get_sprnorm_vecchia()
  R <- 4000
  set.seed(1)
  sim_exact <- sprnorm(spcov_params_val, samples = R, data = exdata, xcoord = xcoord, ycoord = ycoord, local = FALSE)
  set.seed(1)
  sim_vecchia_all <- sprnorm(spcov_params_val, samples = R, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "vecchia", method = "all"))
  expect_equal(rowMeans(sim_exact), rowMeans(sim_vecchia_all), tolerance = 0.1)
  expect_equal(apply(sim_exact, 1, sd), apply(sim_vecchia_all, 1, sd), tolerance = 0.15)

  # neighbor-selection rules and distance/covariance truncation both run
  expect_vector(sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "vecchia", size = 5, method = "distance"))[, 1])
  expect_vector(sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "vecchia", size = 5, method = "covariance"))[, 1])

  # random effects supported via covmatrix() reuse
  exdata_re <- exdata
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  randcov_params_val <- randcov_params(grp = 0.3)
  sim_re <- sprnorm(spcov_params_val, samples = 20, data = exdata_re, xcoord = xcoord, ycoord = ycoord, randcov_params = randcov_params_val, local = list(approximation = "vecchia", size = 10))
  expect_true(all(is.finite(sim_re)))

  # size_base/size_new/reorder_base/kmeans_new (renamed from reorder/kmeans)
  # all still work with an explicit low-rank approximation
  expect_vector(sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "low-rank", size_base = 50, size_new = 20, reorder_base = "random", kmeans_new = FALSE))[, 1])

  # invalid local$approximation errors informatively
  expect_error(sprnorm(spcov_params_val, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(approximation = "bogus")), "local\\$approximation must be")
})

test_that("all spatial simulation families default to Vecchia for local simulation", {
  data <- data.frame(cx = seq(0, 1, length.out = 12), cy = rep(c(0, .3, .1), 4))
  params <- spcov_params("exponential", de = .4, ie = .2, range = .3)
  internal <- data.frame(...xcoord... = data$cx, ...ycoord... = data$cy)
  for (local in list(TRUE, list())) {
    settings <- get_local_list_simulation(local, 12, internal)
    expect_identical(settings$approximation, "vecchia")
    expect_equal(settings$size, 30)
    for (sim in list(sprnorm, sprpois, sprbinom, sprnbinom, sprbeta, sprgamma, sprinvgauss)) {
      set.seed(1)
      default <- sim(params, data = data, xcoord = cx, ycoord = cy, samples = 3, local = local)
      set.seed(1)
      explicit <- sim(params, data = data, xcoord = cx, ycoord = cy, samples = 3,
        local = list(approximation = "vecchia"))
      expect_identical(default, explicit)
    }
  }
  expect_message(settings <- get_local_list_simulation(NULL, 5001, internal), "local = TRUE")
  expect_identical(settings$approximation, "vecchia")
  expect_identical(get_local_list_simulation(FALSE, 12, internal)$method_base, "all")
  expect_error(get_local_list_simulation(list(size_base = 5), 12, internal), "approximation.*low-rank")
})

test_that("Vecchia ordering handles coincident and empty coordinate sets", {
  expect_identical(default_vecchia_ordering(c(0, 0), c(1, 1)), "grts")
  expect_identical(default_vecchia_ordering(c(0, 0), c(1, 2)), "maxmin")
  expect_error(get_decorrelate_order("maxmin", c(0, 0), c(1, 1)), "coincident")
  for (ordering in c("maxmin", "grts", "none", "random")) {
    expect_identical(conditional_vecchia_order(matrix(numeric(), 0, 2), ordering), integer())
    expect_identical(conditional_vecchia_order(matrix(0, 1, 2), ordering), 1L)
  }
})

test_that("Vecchia ordering rejects incorrect orders", {
  skip_if_not_installed("GPvecchia")
  bad_order <- c(1, 1)
  local_mocked_bindings(order_maxmin_exact = function(...) bad_order, .package = "GPvecchia")
  for (bad in list(c(1, 1), c(1, 2.5), c(1, NA), c(1, Inf), 1, c("1", "2"))) {
    bad_order <- bad
    expect_error(get_decorrelate_order("maxmin", c(0, 1), c(0, 1)), "permutation")
  }
})
