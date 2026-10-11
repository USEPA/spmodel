skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

fitted_randcov_fixture <- function(glm = FALSE, areal = FALSE, partition = TRUE, local = FALSE) {
  set.seed(1)
  n <- 48L
  d <- data.frame(cx = runif(n), cy = runif(n), x = rnorm(n), off = seq(-.2, .2, length.out = n),
    g = factor(rep(c("a", "c", "b"), 16), levels = c("c", "b", "a")),
    part = factor(rep(c("A", "B"), each = 24), levels = c("B", "A")))
  d$y <- if (glm) rpois(n, exp(.5 + .2 * d$x + d$off)) else rnorm(n) + d$x + d$off
  args <- list(formula = y ~ x + offset(off), data = d, random = ~ g,
    randcov_initial = randcov_initial(g = .4, known = "given"))
  if (partition) args$partition_factor <- ~ part
  if (glm) args$family <- "poisson"
  if (areal) {
    W <- matrix(0, n, n)
    W[cbind(1:(n - 1L), 2:n)] <- 1
    W <- W + t(W)
    args$W <- W
    args$spcov_initial <- spcov_initial("car", de = .3, ie = .2, range = .25, known = "given")
  } else {
    args$xcoord <- "cx"
    args$ycoord <- "cy"
    args$local <- local
    args$spcov_initial <- spcov_initial("exponential", de = .3, ie = .2, range = .25, known = "given")
  }
  fit <- do.call(if (areal) { if (glm) spgautor else spautor } else { if (glm) spglm else splm }, args)
  list(fit = fit, data = d)
}

test_that("random-effect fitted values are per true partition for all four models", {
  for (glm in c(FALSE, TRUE)) for (areal in c(FALSE, TRUE)) for (partition in c(FALSE, TRUE)) {
    f <- fitted_randcov_fixture(glm, areal, partition)
    fit <- f$fit
    d <- f$data
    pairs <- if (partition) expand.grid(part = levels(d$part), g = levels(d$g), stringsAsFactors = FALSE) else
      data.frame(g = levels(d$g))
    Z <- vapply(seq_len(NROW(pairs)), function(j) as.numeric(d$g == pairs$g[j] &
      if (partition) d$part == pairs$part[j] else TRUE), numeric(NROW(d)))
    labels <- paste0("g", pairs$g)
    if (partition) labels <- paste0(labels, ":part", pairs$part)
    w <- if (glm) fitted(fit, type = "link") else d$y
    r <- w - d$off - as.numeric(model.matrix(fit) %*% coef(fit))
    target <- setNames(.4 * as.numeric(crossprod(Z, solve(as.matrix(covmatrix(fit)), r))), labels)
    actual <- fitted(fit, type = "randcov")[["1 | g"]]
    expect_equal(actual, target, tolerance = if (glm) 1e-5 else 1e-9)
    set.seed(1)
    draws <- conditional(fit, output = "random", samples = 2000)
    expect_identical(rownames(draws), paste0("1 | g:", labels))
    expect_identical(attr(draws, "effect_info")$level, labels)
    expect_lt(max(abs(rowMeans(draws) - target) / (apply(draws, 1, sd) / sqrt(NCOL(draws)))), 5)
  }
})

test_that("local fitted effects average each partition pair across fitting blocks", {
  for (glm in c(FALSE, TRUE)) {
    exact <- fitted_randcov_fixture(glm)
    one <- fitted_randcov_fixture(glm, local = list(index = rep(1L, 48)))
    expect_equal(fitted(one$fit, type = "randcov"), fitted(exact$fit, type = "randcov"),
      tolerance = if (glm) 1e-4 else 1e-9)
    blocks <- rep(rep(1:2, each = 12), 2)
    f <- fitted_randcov_fixture(glm, local = list(index = blocks))
    d <- f$data
    fit <- f$fit
    pairs <- expand.grid(part = levels(d$part), g = levels(d$g), stringsAsFactors = FALSE)
    Z <- vapply(seq_len(NROW(pairs)), function(j) as.numeric(d$g == pairs$g[j] & d$part == pairs$part[j]), numeric(NROW(d)))
    w <- if (glm) fitted(fit, type = "link") else d$y
    r <- w - d$off - as.numeric(model.matrix(fit) %*% coef(fit))
    Sigma <- as.matrix(covmatrix(fit))
    values <- vapply(split(seq_len(NROW(d)), blocks), function(i) {
      .4 * as.numeric(crossprod(Z[i, , drop = FALSE], solve(Sigma[i, i], r[i])))
    }, numeric(NCOL(Z)))
    target <- rowSums(values) / pmax(rowSums(values != 0), 1L)
    expect_equal(unname(fitted(fit, type = "randcov")[[1L]]), target, tolerance = if (glm) 1e-5 else 1e-9)
  }
})

test_that("new effects follow observed factor levels and retain unused fitted zeros", {
  f <- fitted_randcov_fixture(partition = FALSE)
  f$data$g <- factor(f$data$g, levels = c("c", "b", "a", "unused"))
  fit <- update(f$fit, data = f$data)
  expect_identical(names(fitted(fit, type = "randcov")[[1L]]), c("gc", "gb", "ga", "gunused"))
  expect_equal(unname(fitted(fit, type = "randcov")[[1L]]["gunused"]), 0)
  new <- f$data[1:3, ]
  new$g <- c("z", "b", "d")
  draws <- conditional(fit, new, samples = 2, output = "random")
  expect_identical(rownames(draws), paste0("1 | g:g", c("c", "b", "a", "z", "d")))
  expect_identical(vapply(attr(draws, "effect_info")$level_tuple, function(x) as.character(x$g), character(1)),
    c("c", "b", "a", "z", "d"))
})

test_that("random output cannot silently substitute zero draws", {
  f <- fitted_randcov_fixture()
  draw <- get_conditional_random_lm
  local_mocked_bindings(get_conditional_random_lm = function(...) {
    value <- draw(...)
    value$random <- NULL
    value
  })
  expect_error(conditional(f$fit, output = "random", samples = 2), "Requested random-effect draws were not generated.", fixed = TRUE)
})

test_that("crossed slopes use the same interaction labels in fitted and simulated effects", {
  for (glm in c(FALSE, TRUE)) {
    f <- fitted_randcov_fixture(glm, partition = FALSE)
    d <- f$data
    d$h <- factor(rep(c("u", "v"), each = 6, length.out = NROW(d)), levels = c("v", "u"))
    fit <- update(f$fit, data = d, random = ~ (1 | g) + (0 + x | g:h),
      randcov_initial = randcov_initial("1 | g" = .4, "0 + x | g:h" = .2, known = "given"))
    Z <- cbind(model.matrix(~ 0 + g, d), d$x * model.matrix(~ 0 + g:h, d))
    w <- if (glm) fitted(fit, type = "link") else d$y
    r <- w - d$off - as.numeric(model.matrix(fit) %*% coef(fit))
    target <- c(rep(.4, 3), rep(.2, 6)) * as.numeric(crossprod(Z, solve(as.matrix(covmatrix(fit)), r)))
    fitted_effects <- fitted(fit, type = "randcov")
    expect_equal(unname(unlist(fitted_effects)), target, tolerance = if (glm) 1e-5 else 1e-9)
    labels <- unlist(Map(function(term, values) paste0(term, ":", names(values)), names(fitted_effects), fitted_effects), use.names = FALSE)
    expect_identical(rownames(conditional(fit, output = "random", samples = 2)), labels)
  }
})
