skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

test_that("exhaustive simulation neighbors respect order, partitions, and covariance", {
  set.seed(19)
  coordinates <- list(matrix(runif(80), ncol = 2),
    as.matrix(expand.grid(x = 1:6, y = 1:6)), matrix(0, 40, 2),
    cbind(seq_len(40), 0), matrix(c(0, 0), nrow = 1))
  parameters <- lapply(c("exponential", "spherical", "cosine", "wave", "jbessel"),
    function(type) spcov_params(type, de = 1, ie = .2, range = .2))
  for (xy in coordinates) for (params in parameters) for (method in c("distance", "covariance")) {
    n <- nrow(xy)
    partition <- rep(1:3, length.out = n)
    for (size in c(1L, 5L, n)) {
      seed <- .Random.seed
      actual <- sprnorm_vecchia_neighbors(xy, size, method, params, partition)
      expect_identical(.Random.seed, seed)
      distance <- as.matrix(dist(xy))
      for (i in seq_len(n)) {
        pool <- which(seq_len(n) < i & partition == partition[i])
        if (length(pool) > size) {
          score <- if (method == "covariance" && inherits(params, c("cosine", "wave", "jbessel"))) {
            -abs(as.numeric(cov_vector(params, distance[i, pool])))
          } else distance[i, pool]
          pool <- pool[order(score, pool)[seq_len(size)]]
        }
        expect_identical(actual[[i]], pool)
      }
    }
  }
})

test_that("distant sites retain shared random effects with small neighborhoods", {
  d <- data.frame(cx = seq_len(40), cy = 0, g = rep(c("A", "B"), 20),
    part = rep(c("P", "Q"), each = 20), x = seq(-1, 1, length.out = 40))
  params <- spcov_params("exponential", de = .4, ie = .2, range = .001)
  set.seed(231)
  draws <- sprnorm(params, data = d, xcoord = cx, ycoord = cy, samples = 12000,
    randcov_params = c("1 | g" = .5, "x | g" = .3), partition_factor = ~ part,
    local = list(size = 1, ordering = "none"))
  target <- (.5 + .3 * tcrossprod(d$x)) * outer(d$g, d$g, "==") * outer(d$part, d$part, "==") + diag(.6, 40)
  expect_lt(max(abs(cov(t(draws)) - target)), .08)
  expect_lt(max(abs(rowMeans(draws))), .05)
})

test_that("anisotropic full-neighbor simulation reproduces the exact covariance", {
  set.seed(58)
  d <- data.frame(cx = runif(20), cy = runif(20), g = factor(rep(1:4, 5)),
    part = c(rep("A", 9), rep("B", 10), "singleton"))
  params <- spcov_params("exponential", de = .7, ie = .3, range = .3, rotate = .6, scale = .2)
  transformed <- transform_anis(d, "cx", "cy", .6, .2)
  distance <- as.matrix(dist(cbind(transformed$xcoord_val, transformed$ycoord_val)))
  target <- (.7 * exp(-distance / .3) + diag(.3, 20) + .4 * outer(d$g, d$g, "==")) *
    outer(d$part, d$part, "==")
  set.seed(59)
  draws <- sprnorm(params, data = d, xcoord = cx, ycoord = cy, samples = 15000,
    randcov_params = c(g = .4), partition_factor = ~ part,
    local = list(size = 19, ordering = "random"))
  expect_lt(max(abs(cov(t(draws)) - target)), .07)
})

test_that("truncated anisotropic neighborhoods use transformed distances", {
  set.seed(64)
  n <- 20
  d <- data.frame(cx = runif(n), cy = runif(n), response = seq_len(n))
  object <- splm(response ~ 1, data = d, xcoord = cx, ycoord = cy,
    spcov_initial = spcov_initial("exponential", de = .7, ie = .3,
      range = .3, rotate = .6, scale = .2, known = "given"),
    local = list(var_adjust = "none"))
  ord <- sample(n)
  transformed <- transform_anis(d, "cx", "cy", .6, .2)
  distance <- as.matrix(dist(cbind(transformed$xcoord_val, transformed$ycoord_val)[ord, ]))
  raw_distance <- as.matrix(dist(d[ord, c("cx", "cy")]))
  neighbors <- function(distance) lapply(seq_len(n), function(i) {
    pool <- seq_len(i - 1L)
    if (length(pool) > 3L) pool <- pool[order(distance[i, pool], pool)[1:3]]
    pool
  })
  pools <- neighbors(distance)
  expect_false(identical(pools, neighbors(raw_distance)))
  covariance <- as.matrix(covmatrix(object))[ord, ord]
  factor <- matrix(0, n, n)
  for (i in seq_len(n)) {
    pool <- pools[[i]]
    weights <- if (length(pool)) solve(covariance[pool, pool, drop = FALSE], covariance[pool, i]) else numeric()
    factor[i, ] <- as.numeric(crossprod(weights, factor[pool, , drop = FALSE]))
    factor[i, i] <- sqrt(covariance[i, i] - sum(weights * covariance[pool, i]))
  }
  sampler <- get_sprnorm_vecchia
  environment(sampler) <- list2env(list(rnorm = function(count, ...) {
    stopifnot(count == n * n)
    as.numeric(diag(n))
  }),
    parent = environment(get_sprnorm_vecchia))
  actual <- sampler(object, list(order = ord, method = "distance", size = 3), n)
  expect_equal(actual[ord, ], factor, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("full-neighbor simulation has the exact covariance without Monte Carlo error", {
  set.seed(65)
  n <- 24
  d <- data.frame(cx = runif(n), cy = runif(n), response = seq_len(n),
    g = factor(rep(1:4, 6)), h = factor(rep(1:3, 8)), slope = rnorm(n),
    part = c(rep("A", 11), rep("B", 12), "singleton"))
  coincident <- d
  coincident[13:16, c("cx", "cy")] <- coincident[1:4, c("cx", "cy")]
  initial <- function(...) spcov_initial("exponential", de = 1, ie = .2, range = .3, known = "given", ...)
  cases <- list(
    list(data = d, spcov_initial = initial()),
    list(data = d, spcov_initial = initial(rotate = .6, scale = .3)),
    list(data = d, spcov_initial = spcov_initial("matern", de = 1, ie = .2, range = .3, extra = 1.2, known = "given")),
    list(data = d, spcov_initial = spcov_initial("wave", de = 1, ie = .2, range = .1, known = "given")),
    list(data = d, spcov_initial = initial(), random = ~ (1 | g) + (0 + slope | h),
      randcov_initial = randcov_initial(c(`1 | g` = .5, `0 + slope | h` = .3), known = "given"), partition_factor = ~ part),
    list(data = d, spcov_initial = initial(), random = ~ (1 | g) + (0 + slope | g),
      randcov_initial = randcov_initial(c(`1 | g` = .5, `0 + slope | g` = .3), known = "given"), partition_factor = ~ part),
    list(data = coincident, spcov_initial = initial()),
    list(data = coincident, spcov_initial = spcov_initial("exponential", de = 1, ie = 0, range = .3, known = "given"))
  )
  for (case in cases) {
    object <- do.call(splm, c(list(formula = response ~ 1, xcoord = "cx", ycoord = "cy",
      local = list(var_adjust = "none")), case))
    randcov <- coef(object, type = "randcov")
    q <- if (length(randcov)) ncol(get_randcov_design(object$obdata, randcov, object$partition_factor)$observed) else 0L
    sampler <- get_sprnorm_vecchia
    random_sampler <- draw_sprnorm_random
    draw_number <- 0L
    env <- new.env(parent = environment(get_sprnorm_vecchia))
    env$rnorm <- function(count, ...) {
      draw_number <<- draw_number + 1L
      values <- if (draw_number == 1L) cbind(diag(n), matrix(0, n, q)) else cbind(matrix(0, q, n), diag(q))
      stopifnot(length(values) == count)
      as.numeric(values)
    }
    environment(sampler) <- env
    environment(random_sampler) <- env
    env$draw_sprnorm_random <- random_sampler
    draws <- sampler(object, list(order = sample(n), method = "distance", size = n - 1L), n + q)
    expect_equal(tcrossprod(draws), as.matrix(covmatrix(object)), tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("simulation skips unused coefficient covariance adjustments", {
  local_mocked_bindings(get_W_ij = function(...) stop("unused coefficient adjustment"))
  set.seed(66)
  d <- data.frame(cx = runif(240), cy = runif(240))
  params <- spcov_params("exponential", de = .7, ie = .3, range = .3)
  for (local in list(FALSE, list(size = 3, ordering = "none"))) {
    expect_length(sprnorm(params, data = d, xcoord = cx, ycoord = cy, local = local), 240)
  }
})

test_that("iid simulation separates slopes and partition-specific effects without dense algebra", {
  d <- data.frame(x = seq(-1, 1, length.out = 24), g = factor(rep(1:3, 8)), part = factor(rep(1:2, each = 12)))
  local_mocked_bindings(cov_matrix = function(...) stop("dense covariance"),
    partition_matrix = function(...) stop("dense partition"),
    splm = function(...) stop("unnecessary model fit"))
  target <- (.4 + .3 * tcrossprod(d$x)) * outer(d$g, d$g, "==") * outer(d$part, d$part, "==") + diag(.2, 24)
  parameters <- list(spcov_params("none", ie = .2), spcov_params("ie", ie = .2),
    spcov_params("exponential", de = 0, ie = .2, range = .3))
  reference <- NULL
  for (params in parameters) for (local in list(FALSE, TRUE)) {
    set.seed(714)
    draws <- sprnorm(params, data = d, randcov_params = c("1 | g" = .4, "x | g" = .3),
      partition_factor = ~ part, local = local, samples = 15000)
    if (is.null(reference)) reference <- draws else expect_identical(draws, reference)
    expect_lt(max(abs(cov(t(draws)) - target)), .05)
  }
  large <- data.frame(g = factor(seq_len(20000)), part = factor(rep(1:2, 10000)))
  expect_length(sprnorm(parameters[[1]], data = large, randcov_params = c(g = .4),
    partition_factor = ~ part), 20000)
  zero <- spcov_params("none", ie = 0)
  expect_identical(sprnorm(zero, data = d, mean = 2), rep(2, 24))
  expect_true(all(is.finite(sprnorm(zero, data = d, randcov_params = c(g = .4)))))
  bad <- d
  bad$g <- as.numeric(bad$g)
  expect_error(sprnorm(parameters[[1]], data = bad, randcov_params = c(g = .4)), "categorical or factor")
})

test_that("all-neighbor dispatch is exact and avoids ordering", {
  d <- data.frame(cx = rep(1:4, each = 2), cy = 0, g = factor(rep(1:2, 4)))
  params <- spcov_params("exponential", de = .4, ie = .2, range = .3)
  local_mocked_bindings(get_decorrelate_order = function(...) stop("ordering called"))
  for (sim in list(sprnorm, sprpois, sprbinom, sprnbinom, sprbeta, sprgamma, sprinvgauss)) {
    set.seed(6)
    exact <- sim(params, data = d, xcoord = cx, ycoord = cy,
      randcov_params = c(g = .3), samples = 3, local = FALSE)
    set.seed(6)
    all <- sim(params, data = d, xcoord = cx, ycoord = cy,
      randcov_params = c(g = .3), samples = 3, local = list(method = "all"))
    expect_identical(all, exact)
  }
  internal <- setNames(d[1:2], c("...xcoord...", "...ycoord..."))
  expect_error(get_local_list_simulation(list(method = "all", size = -1), 8, internal), "size must")
})
