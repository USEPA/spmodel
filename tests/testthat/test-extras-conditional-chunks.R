skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

chunk_conditional_fixture <- function(glm = FALSE, random = FALSE) {
  set.seed(1)
  d <- data.frame(cx = runif(22), cy = runif(22), x = rnorm(22),
    g = factor(rep(1:2, 11)), off = seq(-.2, .2, length.out = 22))
  d$y <- if (glm) rpois(22, exp(.3 + .2 * d$x + d$off)) else rnorm(22) + d$off
  args <- list(formula = y ~ x + offset(off), data = d[1:18, ],
    xcoord = "cx", ycoord = "cy",
    spcov_initial = spcov_initial("exponential", de = .3, ie = .2, range = .3, known = "given"))
  if (glm) args$family <- "poisson"
  if (random) {
    args$random <- ~ (1 | g)
    args$randcov_initial <- randcov_initial(g = .3, known = "given")
  }
  list(fit = do.call(if (glm) spglm else splm, args), new = d[19:22, ])
}

test_that("local LM and GLM draws bound working columns and reuse preparations", {
  normal_sizes <- integer()
  covariance_calls <- 0L
  prediction_columns <- integer()
  covariance <- covmatrix
  prediction <- draw_conditional_vecchia
  local_mocked_bindings(
    rnorm = function(n, ...) {
      normal_sizes <<- c(normal_sizes, n)
      stats::rnorm(n, ...)
    },
    covmatrix = function(...) {
      covariance_calls <<- covariance_calls + 1L
      covariance(...)
    },
    draw_conditional_vecchia = function(context, base_val, samples) {
      prediction_columns <<- c(prediction_columns, NCOL(base_val))
      prediction(context, base_val, samples)
    })
  for (glm in c(FALSE, TRUE)) for (random in c(FALSE, TRUE)) {
    f <- chunk_conditional_fixture(glm, random)
    for (approximation in c("vecchia", "low-rank")) {
      settings <- if (approximation == "vecchia") list(size = 5, ordering = "none") else
        list(approximation = "low-rank", size_base = 9, size_new = 2,
          reorder_base = "none", reorder_new = "none")
      preparation_counts <- integer()
      for (chunk in c(1L, 3L, 20L)) {
        settings$chunk_size <- chunk
        normal_sizes <- prediction_columns <- integer()
        covariance_calls <- 0L
        output <- c("newdata", "beta", "object")
        if (glm && approximation == "vecchia") output <- c(output, "latent")
        if (random && approximation == "vecchia") output <- c(output, "random")
        set.seed(1)
        ans <- conditional(f$fit, f$new, samples = 7, output = output, local = settings)
        expect_equal(dim(ans$newdata), c(4, 7))
        expect_true(all(is.finite(ans$newdata)))
        expect_equal(dim(ans$beta), c(2, 7))
        expect_lte(max(normal_sizes), f$fit$n * min(chunk, 7))
        if (approximation == "vecchia") {
          expect_equal(prediction_columns, pmin(chunk, 7L - seq.int(0L, 6L, by = chunk)))
        }
        preparation_counts <- c(preparation_counts, covariance_calls)
      }
      expect_equal(preparation_counts, rep(preparation_counts[1], 3))
    }
  }
})

test_that("LM chunks preserve seeded draws, offsets and prediction row order", {
  f <- chunk_conditional_fixture()
  for (settings in list(list(size = 5, ordering = "coordinate"),
      list(approximation = "low-rank", size_base = 9, size_new = 2,
        reorder_base = "none", reorder_new = "random", kmeans_new = FALSE))) {
    settings$chunk_size <- 20
    set.seed(1)
    reference <- conditional(f$fit, f$new, samples = 7, output = "all", local = settings)
    for (chunk in c(1, 3)) {
      settings$chunk_size <- chunk
      set.seed(1)
      ans <- conditional(f$fit, f$new, samples = 7, output = "all", local = settings)
      expect_equal(ans, reference, tolerance = 1e-12)
    }
  }
})

test_that("exact dispatch ignores chunk sizes with and without random effects", {
  effect_columns <- integer()
  effect_draw <- draw_conditional_random
  local_mocked_bindings(draw_conditional_random = function(context, residual) {
    effect_columns <<- c(effect_columns, NCOL(residual))
    effect_draw(context, residual)
  })
  for (glm in c(FALSE, TRUE)) for (random in c(FALSE, TRUE)) {
    f <- chunk_conditional_fixture(glm, random)
    set.seed(1)
    reference <- conditional(f$fit, f$new, samples = 7, local = FALSE)
    for (settings in list(list(method = "all", chunk_size = 1),
        list(approximation = "low-rank", method_base = "all", method_new = "all", chunk_size = 1))) {
      effect_columns <- integer()
      set.seed(1)
      expect_identical(conditional(f$fit, f$new, samples = 7, local = settings), reference)
      if (random) expect_identical(effect_columns, 7L)
    }
  }
})

test_that("local chunks support single draws, single predictions and beta-only output", {
  for (glm in c(FALSE, TRUE)) for (random in c(FALSE, TRUE)) {
    f <- chunk_conditional_fixture(glm, random)
    for (settings in list(list(size = 5, ordering = "none", chunk_size = 1),
        list(approximation = "low-rank", size_base = 9, reorder_base = "none", chunk_size = 1))) {
      set.seed(1)
      expect_equal(dim(conditional(f$fit, f$new[1, ], samples = 1, local = settings)), c(1, 1))
      set.seed(1)
      expect_equal(dim(conditional(f$fit, output = "beta", samples = 3, local = settings)), c(2, 3))
      for (size in list(0, NA, Inf, 1.5, c(1, 2))) {
        settings$chunk_size <- size
        expect_error(conditional(f$fit, f$new, samples = 1, local = settings), "chunk_size must")
      }
    }
  }
})

test_that("chunked low-rank LM blocks support parallel preparation", {
  f <- chunk_conditional_fixture()
  settings <- list(approximation = "low-rank", size_base = 9, size_new = 2,
    reorder_base = "none", reorder_new = "none", chunk_size = 3)
  set.seed(1)
  serial <- conditional(f$fit, f$new, samples = 7, output = "all", local = settings)
  settings$parallel <- TRUE
  settings$ncores <- 2
  set.seed(1)
  parallel <- conditional(f$fit, f$new, samples = 7, output = "all", local = settings)
  expect_identical(parallel$beta, serial$beta)
  expect_equal(dim(parallel$newdata), dim(serial$newdata))
  expect_true(all(is.finite(parallel$newdata)))
})

test_that("local GLM response transformations also respect chunk_size", {
  columns <- integer()
  transform <- invlink_conditional
  local_mocked_bindings(invlink_conditional = function(mu_link, ...) {
    columns <<- c(columns, NCOL(mu_link))
    transform(mu_link, ...)
  })
  for (random in c(FALSE, TRUE)) {
    f <- chunk_conditional_fixture(glm = TRUE, random = random)
    for (settings in list(list(size = 5, ordering = "none", chunk_size = 3),
        list(approximation = "low-rank", size_base = 9, reorder_base = "none", chunk_size = 3))) {
      set.seed(1)
      link <- conditional(f$fit, f$new, samples = 7, local = settings)
      set.seed(1)
      columns <- integer()
      response <- conditional(f$fit, f$new, samples = 7, local = settings, type = "response")
      expect_equal(response, exp(link))
      expect_equal(columns, c(3L, 3L, 1L))
      set.seed(1)
      columns <- integer()
      response <- conditional(f$fit, f$new, samples = 7, local = settings, type = "new")
      expect_true(all(is.finite(response) & response >= 0 & response == floor(response)))
      expect_equal(columns, c(3L, 3L, 1L))
    }
  }
})
