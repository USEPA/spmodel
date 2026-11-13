skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

test_that("one-dimensional simulation warns only for supplied y-coordinates", {
  d <- data.frame(cx = seq(0, 1, length.out = 12), cy = 0)
  for (type in c("triangular", "cosine", "exponential")) {
    params <- spcov_params(type, de = .4, ie = .2, range = .3)
    for (sim in list(sprnorm, sprpois, sprbinom, sprnbinom,
      sprbeta, sprgamma, sprinvgauss)) {
      set.seed(2)
      expect_warning(omitted <- sim(params, data = d, xcoord = cx,
        samples = 3), NA)
      set.seed(2)
      expect_warning(null <- sim(params, data = d, xcoord = cx, ycoord = NULL,
        samples = 3), NA)
      expect_identical(null, omitted)
      set.seed(2)
      warning <- if (type == "exponential") NA else "Ignoring y-coordinate"
      expect_warning(supplied <- sim(params, data = d, xcoord = cx, ycoord = cy,
        samples = 3), warning)
      expect_identical(supplied, omitted)
      expect_true(all(is.finite(omitted)))
    }
  }
})

test_that("sf simulation preserves one-dimensional covariance", {
  set.seed(3)
  data <- data.frame(cx = runif(40), cy = runif(40))
  geometry <- sf::st_as_sf(data, coords = c("cx", "cy"))
  for (type in c("triangular", "cosine")) {
    params <- spcov_params(type, de = 1, ie = .3, range = .3)
    for (sim in list(sprnorm, sprpois)) {
      set.seed(3)
      expect_warning(actual <- sim(params, data = geometry, samples = 3),
        "Ignoring y-coordinate")
      set.seed(3)
      expected <- sim(params, data = data, xcoord = cx, samples = 3)
      expect_equal(actual, expected, tolerance = 1e-12)
    }
  }
})

test_that("shared ordering handles coincident coordinate sets", {
  expect_identical(default_vecchia_ordering(c(0, 0), c(1, 1)), "grts")
  expect_identical(default_vecchia_ordering(c(0, 0), c(1, 2)), "maxmin")
  expect_error(get_decorrelate_order("maxmin", c(0, 0), c(1, 1)), "coincident")
})

test_that("shared ordering rejects incorrect orders", {
  skip_if_not_installed("GPvecchia")
  bad_order <- c(1, 1)
  local_mocked_bindings(order_maxmin_exact = function(...) bad_order, .package = "GPvecchia")
  for (bad in list(c(1, 1), c(1, 2.5), c(1, NA), c(1, Inf), 1, c("1", "2"))) {
    bad_order <- bad
    expect_error(get_decorrelate_order("maxmin", c(0, 1), c(0, 1)), "permutation")
  }
})
