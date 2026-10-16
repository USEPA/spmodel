skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set SPMODEL_RUN_EXTRAS = 'true' to run the extras suite")

test_that("coordinate availability is independent of covariance dimension", {
  d <- data.frame(y = 1:12, cx = 1:12, cy = rep(0:1, 6))
  for (type in c("none", "ie")) {
    params <- spcov_initial(type, ie = .5, known = "given")
    for (coords in list(c("cx", "cy"), list("cx", NULL), list(NULL, "cy"), list(NULL, NULL))) {
      value <- get_point_ref_coords(d, params, "reml", coords[[1]], coords[[2]])
      expect_identical(value$coordinates_supplied, any(lengths(coords) > 0))
      expect_equal(value$dim_coords, 0)
    }
    value <- get_point_ref_coords(sf::st_as_sf(d, coords = c("cx", "cy")), params, "reml", NULL, NULL)
    expect_true(value$coordinates_supplied)
    expect_equal(value$dim_coords, 0)
  }
})

test_that("coordinate-free estimation chooses random groups and respects explicit indices", {
  d <- data.frame(cx = rep(0, 12), cy = 0)
  for (local in list(TRUE, list(), list(size = 4))) {
    set.seed(2)
    actual <- get_local_list_estimation(local, d, "cx", "cy", 12, NULL, FALSE)
    expect_identical(actual$method, "random")
    set.seed(2)
    expect_identical(actual$index, sample(rep(seq_len(actual$groups), times = actual$size)[1:12]))
  }
  expect_error(get_local_list_estimation(list(method = "kmeans"), d, "cx", "cy", 12, NULL, FALSE),
    'kmeans.*requires coordinates.*random')
  index <- rep(1:3, each = 4)
  expect_identical(get_local_list_estimation(list(index = index, method = "kmeans"),
    d, "cx", "cy", 12, NULL, FALSE)$index, index)
  expect_identical(get_local_list_estimation(FALSE, d, "cx", "cy", 12, NULL, FALSE)$index, rep(1, 12))
})

test_that("automatic block thresholds use coordinate availability", {
  for (threshold in c(3000, 5000)) for (n in threshold + c(-1, 0, 1)) {
    d <- data.frame(cx = rep(0, n), cy = 0)
    args <- list(local = NULL, obdata = d, xcoord = "cx", ycoord = "cy", n = n,
      partition_factor = NULL, n_threshold = threshold, local_message = "automatic local",
      X = matrix(1, n, 1), y = matrix(seq_len(n)), coordinates_supplied = FALSE)
    if (n > threshold) {
      expect_message(value <- do.call(build_local_partition_lists, args), "automatic local")
      expect_identical(value$local$method, "random")
      expect_gt(length(unique(value$local$index)), 1)
    } else {
      expect_message(value <- do.call(build_local_partition_lists, args), NA)
      expect_equal(length(unique(value$local$index)), 1)
    }
  }
})

local_coordinate_fixture <- function(n) {
  data.frame(x = sin(seq_len(n)), y = rep(c(1, 3, 2, 4, 2), length.out = n),
    off = .1 * cos(seq_len(n)), g = rep(letters[1:5], length.out = n),
    cx = (seq_len(n) * sqrt(2)) %% 1, cy = (seq_len(n) * sqrt(3)) %% 1)
}

test_that("coordinate-free exported fits work above explicit thresholds", {
  d <- local_coordinate_fixture(3001)
  for (type in c("none", "ie")) {
    params <- spcov_initial(type, ie = .5, known = "given")
    for (local in list(NULL, TRUE)) {
      set.seed(2)
      fit <- suppressMessages(spglm(y ~ x + offset(off), d, family = "poisson",
        spcov_initial = params, local = local))
      expect_true(all(is.finite(coef(fit))))
      expect_equal(fit$n, 3001)
    }
    expect_error(spglm(y ~ x, d, family = "poisson", spcov_initial = params,
      local = list(method = "kmeans")), 'kmeans.*requires coordinates.*random')
  }
  d <- local_coordinate_fixture(5001)
  for (local in list(NULL, TRUE)) {
    set.seed(2)
    fit <- suppressMessages(splm(y ~ x + offset(off), d, random = ~g, ddf = "asymptotic",
      spcov_initial = spcov_initial("none", ie = .5, known = "given"),
      randcov_initial = randcov_initial(g = .3, known = "given"), local = local))
    expect_true(all(is.finite(coef(fit))))
    expect_equal(fit$n, 5001)
  }
})

test_that("supplied coordinates preserve k-means fits", {
  resolve <- get_local_list_estimation
  legacy_default <- function(local, data, xcoord, ycoord, n, partition_factor, coordinates_supplied) {
    expect_true(coordinates_supplied)
    if (isTRUE(local)) local <- list()
    if (is.list(local) && is.null(local$index) && is.null(local$method)) local$method <- "kmeans"
    resolve(local, data, xcoord, ycoord, n, partition_factor)
  }
  for (glm in c(FALSE, TRUE)) for (type in c("none", "ie")) {
    d <- local_coordinate_fixture(if (glm) 3001 else 5001)
    formula <- y ~ x + offset(off)
    random <- ~g
    params <- spcov_initial(type, ie = .5, known = "given")
    randparams <- randcov_initial(g = .3, known = "given")
    for (local in list(NULL, TRUE, list(method = "kmeans"))) {
      run <- function() {
        if (glm) suppressMessages(spglm(formula, d, family = "poisson", spcov_initial = params,
          xcoord = cx, ycoord = cy, local = local)) else
          suppressMessages(splm(formula, d, random = random, randcov_initial = randparams,
            spcov_initial = params, ddf = "asymptotic", xcoord = cx, ycoord = cy, local = local))
      }
      set.seed(2)
      actual <- run()
      set.seed(2)
      reference <- testthat::with_mocked_bindings(run(), get_local_list_estimation = legacy_default,
        .package = "spmodel")
      expect_identical(actual, reference)
    }
  }
})

test_that("stored and supplied coordinate-free prediction rows use the same controls", {
  d <- local_coordinate_fixture(24)
  d$y[21:24] <- NA
  randparams <- randcov_initial(g = .3, known = "given")
  nd <- d[21:24, c("x", "g")]
  for (glm in c(FALSE, TRUE)) for (type in c("none", "ie")) {
    params <- spcov_initial(type, ie = .5, known = "given")
    fit <- if (glm) spglm(y ~ x, d, family = "poisson", random = ~g,
      spcov_initial = params, randcov_initial = randparams) else
      splm(y ~ x, d, random = ~g, spcov_initial = params, randcov_initial = randparams,
        ddf = "asymptotic")
    for (ctl in list(TRUE, list(approximation = "vecchia", ordering = "random"),
      list(approximation = "low-rank", size_base = 10, kmeans_new = FALSE),
      list(approximation = "low-rank", method_base = "all", size_new = 2, kmeans_new = TRUE))) {
      implicit <- tryCatch(conditional(fit, samples = 1, local = ctl), error = conditionMessage)
      expect_match(implicit, "simulation without coordinates")
      expect_error(conditional(fit, nd, samples = 1, local = ctl), implicit, fixed = TRUE)
      expect_error(conditional(fit, fit$newdata, samples = 1, local = ctl), implicit, fixed = TRUE)
    }
    for (ctl in list(FALSE, list(approximation = "vecchia", method = "all"),
      list(approximation = "low-rank", size_base = 10, size_new = 2,
        reorder_base = "random", kmeans_new = FALSE))) {
      set.seed(2)
      implicit <- conditional(fit, samples = 3, local = ctl)
      set.seed(2)
      expect_identical(conditional(fit, nd, samples = 3, local = ctl), implicit)
      expect_identical(dim(implicit), c(4L, 3L))
    }
  }
})

test_that("supplied coordinates remain available for independent random-effect simulation", {
  d <- local_coordinate_fixture(24)
  d$y[21:24] <- NA
  params <- spcov_initial("ie", ie = .5, known = "given")
  randparams <- randcov_initial(g = .3, known = "given")
  for (glm in c(FALSE, TRUE)) for (geometry in c("xy", "sf", "coincident")) {
    data <- d
    if (geometry == "coincident") data$cx <- data$cy <- 0
    args <- list(formula = y ~ x, data = data, random = ~g,
      spcov_initial = params, randcov_initial = randparams)
    if (geometry == "sf") {
      args$data <- sf::st_as_sf(data, coords = c("cx", "cy"))
    } else {
      args$xcoord <- "cx"
      args$ycoord <- "cy"
    }
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fit <- do.call(if (glm) spglm else splm, args)
    nd <- args$data[21:24, ]
    for (ctl in list(TRUE, list(approximation = "low-rank", size_base = 10,
      reorder_base = "grts", kmeans_new = FALSE))) {
      set.seed(2)
      implicit <- conditional(fit, samples = 3, local = ctl)
      set.seed(2)
      expect_equal(conditional(fit, nd, samples = 3, local = ctl), implicit)
      expect_identical(dim(implicit), c(4L, 3L))
    }
    ctl <- list(approximation = "low-rank", method_base = "all", size_new = 2, kmeans_new = TRUE)
    if (geometry == "coincident") {
      expect_error(conditional(fit, samples = 1, local = ctl), "not enough distinct prediction locations.*kmeans_new = FALSE")
    } else {
      expect_identical(dim(conditional(fit, samples = 1, local = ctl)), c(4L, 1L))
    }
  }
})
