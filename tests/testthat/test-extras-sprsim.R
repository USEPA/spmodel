skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

test_that("one-dimensional simulation warns only for supplied y-coordinates", {
  d <- data.frame(cx = seq(0, 1, length.out = 12), cy = 0)
  settings <- list(FALSE, list(size = 4, ordering = "none"))
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



test_that("sprnorm() Vecchia simulation works", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  spcov_params_val <- spcov_params("exponential", de = 1, ie = 0.05, range = 1.5)

  sim1 <- sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(size = 10))
  expect_true(inherits(sim1, "matrix"))
  expect_equal(dim(sim1), c(NROW(exdata), 20))
  expect_true(all(is.finite(sim1)))

  R <- 20
  set.seed(1)
  sim_exact <- sprnorm(spcov_params_val, samples = R, data = exdata, xcoord = xcoord, ycoord = ycoord, local = FALSE)
  set.seed(1)
  sim_vecchia_all <- sprnorm(spcov_params_val, samples = R, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(method = "all"))
  expect_identical(sim_exact, sim_vecchia_all)

  # neighbor-selection rules and distance/covariance truncation both run
  expect_vector(sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(size = 5, method = "distance"))[, 1])
  expect_vector(sprnorm(spcov_params_val, samples = 20, data = exdata, xcoord = xcoord, ycoord = ycoord, local = list(size = 5, method = "covariance"))[, 1])

  # random effects are added to the spatial draws
  exdata_re <- exdata
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  randcov_params_val <- randcov_params(grp = 0.3)
  sim_re <- sprnorm(spcov_params_val, samples = 20, data = exdata_re, xcoord = xcoord, ycoord = ycoord, randcov_params = randcov_params_val, local = list(size = 10))
  expect_true(all(is.finite(sim_re)))


})

test_that("all spatial simulation families default to Vecchia for local simulation", {
  data <- data.frame(cx = seq(0, 1, length.out = 12), cy = rep(c(0, .3, .1), 4))
  params <- spcov_params("exponential", de = .4, ie = .2, range = .3)
  internal <- data.frame(...xcoord... = data$cx, ...ycoord... = data$cy)
  for (local in list(TRUE, list())) {
    settings <- get_local_list_simulation(local, 12, internal)
    expect_false(settings$exact)
    expect_equal(settings$size, 30)
    for (sim in list(sprnorm, sprpois, sprbinom, sprnbinom, sprbeta, sprgamma, sprinvgauss)) {
      set.seed(1)
      default <- sim(params, data = data, xcoord = cx, ycoord = cy, samples = 3, local = local)
      set.seed(1)
      explicit <- sim(params, data = data, xcoord = cx, ycoord = cy, samples = 3,
        local = list())
      expect_identical(default, explicit)
    }
  }
  expect_message(settings <- get_local_list_simulation(NULL, 5001, internal), "local = TRUE")
  expect_false(settings$exact)
  expect_identical(get_local_list_simulation(FALSE, 12, internal)$exact, TRUE)
  expect_error(get_local_list_simulation(list(size_base = 5), 12, internal), "Unknown local argument")
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

test_that("simulation controls reject unknown and malformed local settings", {
  data <- data.frame(...xcoord... = 1:4, ...ycoord... = c(0, 1, 0, 1))
  for (ctl in list(NA, c(TRUE, FALSE), 1, "vecchia", list(3),
      list(size = 3, size = 4), setNames(list(3), NA_character_))) {
    expect_error(get_local_list_simulation(ctl, 4, data), "local must")
  }
  expect_error(get_local_list_simulation(list(unknown = TRUE), 4, data), "Unknown local argument")
  for (size in list(0, -1, NA, Inf, 1.5, c(1, 2), "3")) {
    expect_error(get_local_list_simulation(list(size = size), 4, data), "size must")
  }
})

test_that("sf simulation preserves one-dimensional covariance", {
  set.seed(3)
  data <- data.frame(cx = runif(40), cy = runif(40))
  geometry <- sf::st_as_sf(data, coords = c("cx", "cy"))
  for (type in c("triangular", "cosine")) {
    params <- spcov_params(type, de = 1, ie = .3, range = .3)
    for (sim in list(sprnorm, sprpois)) for (local in list(FALSE, list(size = 10, ordering = "none"))) {
      set.seed(3)
      expect_warning(actual <- sim(params, data = geometry, samples = 3, local = local),
        "Ignoring y-coordinate")
      set.seed(3)
      expected <- sim(params, data = data, xcoord = cx, samples = 3, local = local)
      expect_equal(actual, expected, tolerance = 1e-12)
    }
  }
})
