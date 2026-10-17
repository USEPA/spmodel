skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set SPMODEL_RUN_EXTRAS = 'true' to run the extras suite")

local_coordinate_fixture <- function(n) {
  data.frame(x = sin(seq_len(n)), y = rep(c(1, 3, 2, 4, 2), length.out = n),
    off = .1 * cos(seq_len(n)), g = rep(letters[1:5], length.out = n),
    cx = (seq_len(n) * sqrt(2)) %% 1, cy = (seq_len(n) * sqrt(3)) %% 1)
}

test_that("explicit NULL coordinates match omitted coordinates in independent fits", {
  d <- local_coordinate_fixture(24)
  d$y[21:24] <- NA
  nd <- d[21:24, ]
  for (glm in c(FALSE, TRUE)) for (type in c("none", "ie")) for (random in c(FALSE, TRUE)) {
    args <- list(formula = y ~ x + offset(off), data = d,
      spcov_initial = spcov_initial(type, ie = .5, known = "given"),
      local = list(size = 8, method = "random"))
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    if (random) {
      args$random <- ~g
      args$randcov_initial <- randcov_initial(g = .3, known = "given")
    }
    fun <- if (glm) spglm else splm
    set.seed(2)
    reference <- do.call(fun, args)
    for (coords in list(list(xcoord = NULL), list(ycoord = NULL),
      list(xcoord = NULL, ycoord = NULL))) {
      set.seed(2)
      actual <- do.call(fun, c(args, coords))
      expect_identical(actual[setdiff(names(actual), "call")],
        reference[setdiff(names(reference), "call")])
      expect_false(conditional_has_coordinates(actual, actual$newdata))
      expect_identical(predict(actual, nd, se.fit = TRUE, local = FALSE),
        predict(reference, nd, se.fit = TRUE, local = FALSE))
      for (ctl in list(FALSE, list(approximation = "low-rank", size_base = 10,
        size_new = 2, reorder_base = "random", kmeans_new = FALSE))) {
        set.seed(2)
        expected <- conditional(reference, samples = 3, local = ctl)
        set.seed(2)
        expect_identical(conditional(actual, samples = 3, local = ctl), expected)
        set.seed(2)
        expect_identical(conditional(actual, nd, samples = 3, local = ctl), expected)
      }
      if (random) {
        for (newdata in list(nd, actual$newdata)) {
          expect_error(conditional(actual, newdata, samples = 1, local = TRUE),
            "simulation without coordinates")
        }
        expect_error(conditional(actual, samples = 1, local = TRUE),
          "simulation without coordinates")
      }
    }
  }
})

test_that("explicit NULL coordinates preserve spatial and sf fitting behavior", {
  d <- local_coordinate_fixture(24)
  for (glm in c(FALSE, TRUE)) {
    args <- list(formula = y ~ x, data = d,
      spcov_initial = spcov_initial("exponential", de = .3, ie = .5, range = .4,
        known = "given"))
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fun <- if (glm) spglm else splm
    expect_error(do.call(fun, c(args, list(xcoord = NULL))), "xcoord argument must be specified")
    expect_error(do.call(fun, c(args, list(xcoord = NULL, ycoord = "cy"))),
      "xcoord argument must be specified")
    reference <- do.call(fun, c(args, list(xcoord = "cx")))
    actual <- do.call(fun, c(args, list(xcoord = "cx", ycoord = NULL)))
    expect_identical(actual[setdiff(names(actual), "call")],
      reference[setdiff(names(reference), "call")])
    args$data <- sf::st_as_sf(d, coords = c("cx", "cy"))
    reference <- do.call(fun, args)
    expect_warning(actual <- do.call(fun, c(args, list(xcoord = NULL, ycoord = NULL))), NA)
    expect_identical(actual[setdiff(names(actual), "call")],
      reference[setdiff(names(reference), "call")])
    expect_true(conditional_has_coordinates(actual, args$data))
  }
})

test_that("explicit NULL coordinates preserve unconditional simulation behavior", {
  d <- local_coordinate_fixture(24)
  params <- spcov_params("exponential", de = .3, ie = .5, range = .4)
  for (ctl in list(FALSE, list(approximation = "vecchia", size = 5, ordering = "none"),
    list(approximation = "low-rank", size_base = 10, size_new = 4,
      reorder_base = "none", kmeans_new = FALSE))) {
    set.seed(2)
    reference <- sprnorm(params, data = d, xcoord = cx, samples = 3, local = ctl)
    set.seed(2)
    expect_identical(sprnorm(params, data = d, xcoord = cx, ycoord = NULL,
      samples = 3, local = ctl), reference)
  }
  expect_error(sprnorm(params, data = d, xcoord = NULL), "xcoord argument must be specified")
  for (type in c("none", "ie")) {
    params <- spcov_params(type, ie = .5)
    set.seed(2)
    reference <- sprnorm(params, data = d, samples = 3)
    set.seed(2)
    expect_identical(sprnorm(params, data = d, samples = 3,
      xcoord = NULL, ycoord = NULL), reference)
  }
})

test_that("explicit NULL coordinates preserve empirical diagnostics", {
  d <- local_coordinate_fixture(24)
  without_call <- function(x) {
    attr(x, "call") <- NULL
    x
  }
  for (fun in list(esv, eacf)) {
    args <- list(formula = y ~ x, data = d)
    reference <- do.call(fun, c(args, list(xcoord = "cx")))
    expect_identical(without_call(do.call(fun, c(args, list(xcoord = "cx", ycoord = NULL)))),
      without_call(reference))
    expect_error(do.call(fun, c(args, list(xcoord = NULL))), "xcoord argument must be specified")
    expect_error(do.call(fun, c(args, list(xcoord = "absent", ycoord = NULL))),
      "xcoord argument must match")
    args$dist_matrix <- as.matrix(dist(d["cx"]))
    expect_identical(without_call(do.call(fun, c(args, list(xcoord = NULL, ycoord = NULL)))),
      without_call(do.call(fun, args)))
    args$dist_matrix <- NULL
    args$data <- sf::st_as_sf(d, coords = c("cx", "cy"))
    expect_identical(without_call(do.call(fun, c(args, list(xcoord = NULL, ycoord = NULL)))),
      without_call(do.call(fun, args)))
  }
})

test_that("explicit NULL coordinates preserve decorrelation transforms", {
  d <- local_coordinate_fixture(24)
  for (type in c("none", "ie", "exponential")) for (local in list(FALSE, list(size = 5))) {
    params <- if (type == "exponential") spcov_params(type, de = .3, ie = .5, range = .4) else
      spcov_params(type, ie = .5)
    args <- list(formula = y ~ x, data = d, spcov_params = params,
      ordering = "none", local = local)
    if (type == "exponential") args$xcoord <- "cx"
    reference <- do.call(decorrelate_data, args)
    coords <- if (type == "exponential") list(ycoord = NULL) else list(xcoord = NULL, ycoord = NULL)
    expect_identical(do.call(decorrelate_data, c(args, coords)), reference)
    if (type == "exponential") {
      args$xcoord <- NULL
      expect_error(do.call(decorrelate_data, c(args, list(xcoord = NULL))),
        "xcoord argument must be specified")
    }
    args$xcoord <- NULL
    args$data <- sf::st_as_sf(d, coords = c("cx", "cy"))
    expect_identical(do.call(decorrelate_data, c(args, list(xcoord = NULL, ycoord = NULL))),
      do.call(decorrelate_data, args))
  }
})

test_that("explicit NULL coordinates and one-dimensional decorrelation grids agree", {
  d <- local_coordinate_fixture(24)
  d$zero <- 0
  for (type in c("exponential", "triangular", "cosine", "none", "ie")) {
    args <- list(formula = y ~ x, data = d, spcov_type = type, dense_grid = FALSE)
    reference <- do.call(decorrelate_grid, c(args, list(xcoord = "cx", ycoord = "zero")))
    expect_identical(do.call(decorrelate_grid, c(args, list(xcoord = "cx"))), reference)
    expect_identical(do.call(decorrelate_grid, c(args, list(xcoord = "cx", ycoord = NULL))), reference)
    if (type %in% c("none", "ie")) {
      expect_identical(do.call(decorrelate_grid, args), reference)
      expect_identical(do.call(decorrelate_grid, c(args, list(xcoord = NULL, ycoord = NULL))), reference)
    } else {
      expect_error(do.call(decorrelate_grid, args), "xcoord argument must be specified")
      expect_error(do.call(decorrelate_grid, c(args, list(xcoord = NULL))),
        "xcoord argument must be specified")
    }
    expect_error(do.call(decorrelate_grid, c(args, list(xcoord = "absent"))), "xcoord argument must match")
    expect_error(do.call(decorrelate_grid, c(args, list(xcoord = "cx", ycoord = "absent"))),
      "ycoord argument must match")
    args$data <- sf::st_as_sf(d, coords = c("cx", "cy"))
    expect_identical(do.call(decorrelate_grid, c(args, list(xcoord = NULL, ycoord = NULL))),
      do.call(decorrelate_grid, args))
  }
})

test_that("explicit NULL coordinates preserve decorrelation fits and predictions", {
  skip_if_not_installed("ranger")
  d <- local_coordinate_fixture(24)
  for (type in c("exponential", "none")) for (search in c(FALSE, TRUE)) {
    args <- list(formula = y ~ x, data = d, ordering = "none", local = FALSE,
      num.trees = 10, num.threads = 1)
    if (type == "exponential") args$xcoord <- "cx"
    if (search) {
      args$spcov_type <- type
      args$dense_grid <- FALSE
    } else {
      args$spcov_params <- if (type == "exponential") spcov_params(type, de = .3, ie = .5, range = .4) else
        spcov_params(type, ie = .5)
    }
    set.seed(2)
    reference <- do.call(decorrelate, args)
    coords <- if (type == "exponential") list(ycoord = NULL) else list(xcoord = NULL, ycoord = NULL)
    set.seed(2)
    actual <- do.call(decorrelate, c(args, coords))
    expect_identical(actual[setdiff(names(actual), "call")], reference[setdiff(names(reference), "call")])
    expect_identical(predict(actual, d[1:4, ]), predict(reference, d[1:4, ]))
  }
})

test_that("one-dimensional conditional simulation matches constant-y spatial fits", {
  for (glm in c(FALSE, TRUE)) for (xname in c("cx", ".ycoord")) for (null_y in c(FALSE, TRUE)) {
    d <- local_coordinate_fixture(30)
    d$cy <- 0
    d[[xname]] <- d$cx
    d$y[25:30] <- NA
    nd <- d[25:30, ]
    args <- list(formula = y ~ x + offset(off), data = d, xcoord = xname,
      spcov_initial = spcov_initial("exponential", de = .3, ie = .5, range = .4, known = "given"))
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fun <- if (glm) spglm else splm
    fit <- do.call(fun, c(args, if (null_y) list(ycoord = NULL)))
    reference <- do.call(fun, c(args, list(ycoord = "cy")))
    expect_true(conditional_has_coordinates(fit, nd))
    expect_false(conditional_has_coordinates(fit, nd[setdiff(names(nd), xname)]))
    for (ctl in list(FALSE, TRUE, list(approximation = "vecchia", size = 7, ordering = "none"),
      list(approximation = "low-rank", size_base = 15, size_new = 3,
        reorder_base = "none", kmeans_new = TRUE),
      list(approximation = "low-rank", size_base = 15, size_new = 3,
        reorder_base = "none", kmeans_new = FALSE))) {
      set.seed(2)
      expected <- conditional(reference, nd, samples = 3, local = ctl)
      set.seed(2)
      actual <- conditional(fit, nd, samples = 3, local = ctl)
      expect_equal(actual, expected, tolerance = 1e-12)
      expect_identical(dim(actual), c(6L, 3L))
      set.seed(2)
      expect_identical(conditional(fit, samples = 3, local = ctl), actual)
    }
  }
})

test_that("one-dimensional conditional simulation ignores unused y coordinates", {
  d <- local_coordinate_fixture(24)
  nd <- d[1:6, ]
  nd$cx <- nd$cx + .01
  for (type in c("triangular", "cosine")) {
    params <- spcov_initial(type, de = .3, ie = .5, range = .4, known = "given")
    fit <- splm(y ~ x, d, xcoord = cx, spcov_initial = params, ddf = "asymptotic")
    expect_warning(with_y <- splm(y ~ x, d, xcoord = cx, ycoord = cy,
      spcov_initial = params, ddf = "asymptotic"), "Ignoring y-coordinate")
    for (ctl in list(FALSE, TRUE, list(approximation = "low-rank", size_base = 15,
      size_new = 3, reorder_base = "none", kmeans_new = TRUE))) {
      set.seed(2)
      expected <- conditional(fit, nd, samples = 3, local = ctl)
      set.seed(2)
      expect_identical(conditional(with_y, nd, samples = 3, local = ctl), expected)
      polluted <- nd
      polluted[[fit$ycoord]] <- seq_len(nrow(nd))
      set.seed(2)
      expect_identical(conditional(fit, polluted, samples = 3, local = ctl), expected)
    }
  }
})

test_that("one-dimensional conditional coordinates support geometry and random effects", {
  d <- local_coordinate_fixture(24)
  d$cy <- 0
  nd <- d[1:6, ]
  nd$cx <- nd$cx + .01
  params <- spcov_initial("exponential", de = .3, ie = .5, range = .4, known = "given")
  for (glm in c(FALSE, TRUE)) {
    args <- list(formula = y ~ x + offset(off), data = d, xcoord = "cx", random = ~g,
      spcov_initial = params, randcov_initial = randcov_initial(g = .2, known = "given"))
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fit <- do.call(if (glm) spglm else splm, args)
    for (ctl in list(FALSE, TRUE, list(approximation = "low-rank", size_base = 15,
      size_new = 3, reorder_base = "none", kmeans_new = TRUE))) {
      set.seed(2)
      expected <- conditional(fit, nd, samples = 3, local = ctl)
      set.seed(2)
      expect_equal(conditional(fit, sf::st_as_sf(nd, coords = c("cx", "cy")),
        samples = 3, local = ctl), expected, tolerance = 1e-12)
    }
    expect_error(conditional(fit, nd[setdiff(names(nd), "cx")], local = TRUE), "cx")
  }
})

test_that("one-dimensional conditional coordinates preserve coordinate predictors", {
  d <- local_coordinate_fixture(24)
  d$.conditional_ycoord <- cos(seq_len(nrow(d)))
  nd <- d[1:6, ]
  nd$cx <- nd$cx + .01
  for (glm in c(FALSE, TRUE)) for (type in c("triangular", "cosine")) {
    args <- list(formula = y ~ x + cy + .conditional_ycoord + offset(off), data = d, xcoord = "cx",
      spcov_initial = spcov_initial(type, de = .3, ie = .5, range = .4, known = "given"))
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fun <- if (glm) spglm else splm
    reference <- do.call(fun, args)
    expect_warning(fit <- do.call(fun, c(args, list(ycoord = "cy"))), "Ignoring y-coordinate")
    observed <- fit$obdata
    for (ctl in list(FALSE, TRUE, list(approximation = "low-rank", size_base = 15,
      size_new = 3, reorder_base = "none", kmeans_new = TRUE))) {
      set.seed(2)
      expected <- conditional(reference, nd, samples = 3, local = ctl)
      set.seed(2)
      expect_equal(conditional(fit, nd, samples = 3, local = ctl), expected, tolerance = 1e-12)
    }
    expect_identical(fit$obdata, observed)
  }
})

test_that("one-dimensional covariance fits retain ignored coordinates as data", {
  d <- local_coordinate_fixture(24)
  d$y[21:24] <- NA
  for (glm in c(FALSE, TRUE)) for (type in c("triangular", "cosine")) {
    args <- list(formula = y ~ x + offset(off), data = d, xcoord = "cx",
      spcov_initial = spcov_initial(type, de = .3, ie = .5, range = .4, known = "given"))
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fun <- if (glm) spglm else splm
    reference <- do.call(fun, args)
    expect_warning(fit <- do.call(fun, c(args, list(ycoord = "cy"))), "Ignoring y-coordinate")
    expect_identical(fit$ycoord, reference$ycoord)
    expect_equal(fit$obdata[[fit$ycoord]], rep(0, 20))
    expect_equal(fit$newdata[[fit$ycoord]], rep(0, 4))
    expect_identical(fit$obdata$cy, d$cy[1:20])
    expect_identical(fit$newdata$cy, d$cy[21:24])
    expect_identical(anyDuplicated(names(fit$obdata)), 0L)
    expect_identical(anyDuplicated(names(fit$newdata)), 0L)
    expect_equal(coef(fit), coef(reference))
    expect_equal(coef(fit, type = "spcov"), coef(reference, type = "spcov"))
    expect_equal(logLik(fit), logLik(reference))
    expect_equal(model.matrix(fit), model.matrix(reference))
    expect_equal(fitted(fit), fitted(reference))
    distance <- abs(outer(d$cx[1:20], d$cx[1:20], "-"))
    correlation <- if (type == "triangular") pmax(1 - distance / .4, 0) else cos(distance / .4)
    expect_equal(unname(as.matrix(covmatrix(fit))), .3 * correlation + diag(.5, 20))
    expect_equal(covmatrix(fit), covmatrix(reference))
    for (cov_type in c("obs.pred", "pred.pred")) {
      expect_equal(covmatrix(fit, d[21:24, ], cov_type = cov_type),
        covmatrix(reference, d[21:24, ], cov_type = cov_type))
      expect_equal(covmatrix(fit, fit$newdata, cov_type = cov_type),
        covmatrix(reference, reference$newdata, cov_type = cov_type))
    }
  }
})

test_that("prediction covariance preserves one-dimensional coordinates in newdata", {
  d <- local_coordinate_fixture(24)
  d$part <- factor(rep(1:2, each = 12))
  nd <- d[seq(1, 24, by = 3), ]
  nd$cx <- nd$cx + .01
  for (glm in c(FALSE, TRUE)) for (type in c("triangular", "cosine", "exponential")) {
    args <- list(formula = y ~ x + cy + offset(off), data = d, xcoord = "cx",
      spcov_initial = spcov_initial(type, de = .3, ie = .5, range = .4, known = "given"),
      random = ~g, randcov_initial = randcov_initial(g = .2, known = "given"),
      partition_factor = ~part)
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fun <- if (glm) spglm else splm
    fit <- do.call(fun, args)
    expected <- covmatrix(fit, nd, cov_type = "pred.pred")
    original <- fit$obdata
    collision <- nd
    collision$.ycoord <- nd$cy
    expect_equal(covmatrix(fit, collision, cov_type = "pred.pred"), expected)
    expect_identical(fit$obdata, original)
    expect_identical(collision$.ycoord, nd$cy)
    expect_gt(min(eigen(as.matrix(expected), symmetric = TRUE, only.values = TRUE)$values), 0)

    if (type != "exponential") {
      args$data <- sf::st_as_sf(d, coords = c("cx", "cy"), remove = FALSE)
      args$xcoord <- NULL
      expect_warning(sf_fit <- do.call(fun, args), NA)
      sf_nd <- sf::st_as_sf(nd, coords = c("cx", "cy"), remove = FALSE)
      original <- sf_nd
      expect_equal(unname(as.matrix(covmatrix(sf_fit, sf_nd, cov_type = "pred.pred"))),
        unname(as.matrix(expected)))
      expect_equal(unname(as.matrix(covmatrix(sf_fit, sf_nd, cov_type = "obs.pred"))),
        unname(as.matrix(covmatrix(fit, nd, cov_type = "obs.pred"))))
      expect_identical(sf_nd, original)
    }

    args$data <- d
    names(args$data)[names(args$data) == "cx"] <- ".ycoord"
    args$xcoord <- ".ycoord"
    renamed_fit <- do.call(fun, args)
    renamed <- nd
    names(renamed)[names(renamed) == "cx"] <- ".ycoord"
    renamed$.ycoord2 <- nd$cy
    expect_equal(covmatrix(renamed_fit, renamed, cov_type = "pred.pred"), expected)
    expect_identical(renamed$.ycoord, nd$cx)
    expect_identical(renamed$.ycoord2, nd$cy)
  }
})

test_that("one-dimensional covariance fits agree in prediction and cross validation", {
  d <- local_coordinate_fixture(24)
  d$y[21:24] <- NA
  nd <- d[21:24, setdiff(names(d), "cy")]
  for (glm in c(FALSE, TRUE)) for (type in c("triangular", "cosine")) {
    args <- list(formula = y ~ x + offset(off), data = d, xcoord = "cx",
      spcov_initial = spcov_initial(type, de = .3, ie = .5, range = .4, known = "given"))
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fun <- if (glm) spglm else splm
    reference <- do.call(fun, args)
    expect_warning(fit <- do.call(fun, c(args, list(ycoord = "cy"))), "Ignoring y-coordinate")
    for (ctl in list(FALSE, list(method = "covariance", size = 8),
      list(method = "distance", size = 8))) {
      expect_warning(actual <- predict(fit, nd, se.fit = TRUE,
        interval = "prediction", local = ctl), NA)
      expect_equal(actual, predict(reference, nd, se.fit = TRUE,
        interval = "prediction", local = ctl))
      expect_true(all(is.finite(actual$se.fit)))
      stored <- predict(fit, se.fit = TRUE, interval = "prediction", local = ctl)
      expect_equal(stored, predict(reference, se.fit = TRUE, interval = "prediction", local = ctl))
      expect_equal(lapply(stored, unname), lapply(actual, unname))
      expect_equal(augment(fit, newdata = nd, se_fit = TRUE, local = ctl),
        augment(reference, newdata = nd, se_fit = TRUE, local = ctl))
      expect_equal(loocv(fit, cv_predict = TRUE, se.fit = TRUE, local = ctl),
        loocv(reference, cv_predict = TRUE, se.fit = TRUE, local = ctl))
      set.seed(2)
      expect_warning(expected <- kcv(reference, folds_index = rep(1:4, 5), cv_predict = TRUE,
        se.fit = TRUE, local = ctl), NA)
      set.seed(2)
      expect_warning(actual <- kcv(fit, folds_index = rep(1:4, 5), cv_predict = TRUE,
        se.fit = TRUE, local = ctl), NA)
      expect_equal(actual, expected)
    }
    for (ctl in list(FALSE, TRUE, list(approximation = "low-rank", size_base = 12,
      size_new = 2, reorder_base = "none", kmeans_new = FALSE))) {
      set.seed(2)
      expected <- conditional(reference, nd, samples = 3, local = ctl)
      set.seed(2)
      expect_equal(conditional(fit, nd, samples = 3, local = ctl), expected)
      set.seed(2)
      expect_equal(conditional(fit, samples = 3, local = ctl), expected)
    }
  }
})

test_that("one-dimensional covariance fits preserve coordinate predictors and dot formulas", {
  d <- local_coordinate_fixture(24)
  nd <- d[1:4, ]
  nd$cx <- nd$cx + .01
  for (glm in c(FALSE, TRUE)) for (type in c("triangular", "cosine")) {
    args <- list(formula = y ~ x + cy + offset(off), data = d, xcoord = "cx",
      spcov_initial = spcov_initial(type, de = .3, ie = .5, range = .4, known = "given"))
    if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
    fun <- if (glm) spglm else splm
    reference <- do.call(fun, args)
    expect_warning(fit <- do.call(fun, c(args, list(ycoord = "cy"))), "Ignoring y-coordinate")
    expect_identical(fit$obdata$cy, d$cy)
    expect_equal(model.matrix(fit), model.matrix(reference))
    expect_equal(predict(fit, nd, se.fit = TRUE, local = FALSE),
      predict(reference, nd, se.fit = TRUE, local = FALSE))
    args$data <- d[c("y", "x", "cx", "cy")]
    args$formula <- y ~ x
    reference <- do.call(fun, args)
    args$formula <- y ~ .
    expect_warning(fit <- do.call(fun, c(args, list(ycoord = "cy"))), "Ignoring y-coordinate")
    expect_identical(colnames(model.matrix(fit)), c("(Intercept)", "x"))
    expect_equal(coef(fit), coef(reference))
    expect_equal(predict(fit, nd, se.fit = TRUE, local = FALSE),
      predict(reference, nd, se.fit = TRUE, local = FALSE))
  }
})

test_that("local kcv omits unused coordinates without changing fold results", {
  d <- local_coordinate_fixture(24)
  for (glm in c(FALSE, TRUE)) for (type in c("triangular", "cosine", "exponential", "none", "ie")) {
    for (supplied in c(FALSE, TRUE)) {
      spatial <- !type %in% c("none", "ie")
      params <- if (spatial) spcov_initial(type, de = .3, ie = .5, range = .4, known = "given") else
        spcov_initial(type, ie = .5, known = "given")
      args <- list(formula = y ~ x + offset(off), data = d, spcov_initial = params,
        local = list(index = rep(1:3, each = 8)))
      if (type %in% c("triangular", "cosine")) args$formula <- y ~ x + cy + offset(off)
      if (!spatial) {
        args$random <- ~g
        args$randcov_initial <- randcov_initial(g = .2, known = "given")
      }
      if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
      if (spatial || supplied) args$xcoord <- "cx"
      if (supplied) args$ycoord <- "cy"
      fun <- if (glm) spglm else splm
      if (supplied && type %in% c("triangular", "cosine")) {
        expect_warning(fit <- do.call(fun, args), "Ignoring y-coordinate")
      } else {
        expect_warning(fit <- do.call(fun, args), NA)
      }
      legacy_refit <- function(...) {
        refit_args <- list(...)
        refit_args$ycoord <- fit$ycoord
        do.call(fun, refit_args)
      }
      run <- function() kcv(fit, folds_index = rep(1:4, 6), local = TRUE,
        cv_predict = TRUE, se.fit = TRUE)
      set.seed(2)
      expected <- if (glm) {
        with_mocked_bindings(suppressWarnings(run()), spglm = legacy_refit, .package = "spmodel")
      } else {
        with_mocked_bindings(suppressWarnings(run()), splm = legacy_refit, .package = "spmodel")
      }
      set.seed(2)
      expect_warning(actual <- run(), NA)
      expect_identical(actual, expected)
      expect_true(all(is.finite(actual$cv_predict)) && all(is.finite(actual$se.fit)))
    }
  }
})

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
