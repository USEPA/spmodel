skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

random_conditional_fixture <- function(glm = FALSE, partition = FALSE) {
  set.seed(1)
  d <- data.frame(cx = runif(36), cy = runif(36), x = rnorm(36),
    g = factor(rep(1:4, 9)), h = factor(rep(1:3, each = 12)),
    part = factor(rep(1:2, each = 18)), off = seq(-.2, .2, length.out = 36))
  d$y <- if (glm) rpois(36, exp(.3 + .2 * d$x + d$off)) else rnorm(36) + d$off
  args <- list(formula = y ~ x + offset(off), data = d[1:30, ],
    xcoord = "cx", ycoord = "cy", random = ~ (1 | g),
    spcov_initial = spcov_initial("exponential", de = .3, ie = .2, range = .3, known = "given"),
    randcov_initial = randcov_initial(g = .3, known = "given"))
  if (partition) args$partition_factor <- ~ part
  if (glm) args$family <- "poisson"
  list(fit = do.call(if (glm) spglm else splm, args), new = d[31:36, ], data = d)
}

test_that("spatial prediction covariance retains new character grouping levels", {
  for (partition in c(FALSE, TRUE)) {
    f <- random_conditional_fixture(partition = partition)
    new <- f$new
    new$g <- as.character(new$g)
    new$g[1:4] <- "unseen"
    new$part <- as.character(new$part)
    new$part[1:2] <- "new partition"
    distance <- unname(as.matrix(dist(new[c("cx", "cy")])))
    spatial <- .3 * exp(-distance / .3)
    expected <- spatial + .2 * diag(NROW(new)) + .3 * outer(new$g, new$g, "==")
    if (partition) expected <- expected * outer(new$part, new$part, "==")
    actual <- unname(covmatrix(f$fit, new, cov_type = "pred.pred"))
    expect_equal(actual, expected, tolerance = 1e-12)
    expect_equal(diag(actual), rep(.8, NROW(new)))
    expect_equal(actual[1, 2], spatial[1, 2] + .3)
    if (partition) expect_equal(actual[1, 3], 0)
  }
})

test_that("effect designs preserve covariance, partitions and new levels", {
  for (partition in c(FALSE, TRUE)) {
    fixture <- random_conditional_fixture(partition = partition)
    fit <- fixture$fit
    new <- fixture$new
    new$g <- as.character(new$g)
    new$g[1:2] <- "unseen"
    design <- get_conditional_random_design(fit, new)
    spatial <- conditional_spatial_object(fit)
    expect_equal(as.matrix(covmatrix(fit) - covmatrix(spatial)),
      as.matrix(design$observed %*% Matrix::Diagonal(x = design$variance) %*% t(design$observed)),
      ignore_attr = TRUE)
    expect_equal(sum(!design$info$conditioned), 1)
    expect_identical(design$newdata[1, ], design$newdata[2, ])
    expect_false(anyDuplicated(rownames(design$info)) > 0)
    expect_true(all(c("level_tuple", "partition_tuple") %in% names(design$info)))
    set.seed(1)
    out <- conditional(fit, new, samples = 2000, output = c("random", "newdata"))
    a <- out$random[!design$info$conditioned, ]
    expect_equal(var(as.numeric(a)), .3, tolerance = .03)
    expect_true(all(is.finite(out$newdata)))
    prepared <- prepare_conditional_random(fit, design)
    expect_equal(tcrossprod(prepared$factor), prepared$covariance)
  }
})

test_that("the effect precision agrees with the dense conditional", {
  f <- random_conditional_fixture()
  design <- get_conditional_random_design(f$fit, f$new)
  context <- prepare_conditional_random(f$fit, design)
  Z <- as.matrix(design$observed)
  Gamma <- diag(design$variance)
  Sigma <- as.matrix(covmatrix(f$fit))
  T <- Gamma %*% t(Z) %*% solve(Sigma)
  H <- Gamma - T %*% Z %*% Gamma
  expect_equal(context$covariance, H, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(context$covariance %*% t(context$blocks[[1]]$weights), T,
    ignore_attr = TRUE, tolerance = 1e-10)
  settings <- get_local_list_conditional(list(size = 5, size_random = 10), f$fit, f$new)
  set.seed(1)
  blocks <- prepare_conditional_random(f$fit, design, settings)
  precision <- solve(Gamma)
  S <- as.matrix(covmatrix(conditional_spatial_object(f$fit)))
  for (block in blocks$blocks) {
    j <- block$rows
    precision <- precision + t(Z[j, ]) %*% solve(S[j, j], Z[j, ])
  }
  expect_equal(blocks$covariance, solve(precision), ignore_attr = TRUE)
  expect_equal(sort(unlist(lapply(blocks$blocks, `[[`, "rows"))), seq_len(f$fit$n))
  expect_equal(lengths(conditional_random_blocks(f$fit, 1)), rep(1L, f$fit$n))
  repeated <- f$fit
  repeated$obdata$cx <- repeated$obdata$cy <- 0
  expect_error(conditional_random_blocks(repeated, 10), "Too few distinct")
})

test_that("effect grouping controls follow local fitting conventions", {
  f <- random_conditional_fixture()
  n <- f$fit$n
  design <- get_conditional_random_design(f$fit, f$new)
  Z <- as.matrix(design$observed)
  S <- as.matrix(covmatrix(conditional_spatial_object(f$fit)))
  for (method in c("kmeans", "random")) {
    set.seed(1)
    local <- get_local_list_conditional(list(groups_random = 7, method_random = method), f$fit, f$new)
    expect_equal(local$groups_random, 7)
    expect_equal(local$size_random, ceiling(n / 7))
    expect_length(local$random_blocks, 7)
    expect_equal(sort(unlist(local$random_blocks)), seq_len(n))
    set.seed(1)
    labels <- get_local_estimation_index(list(groups = 7, size = ceiling(n / 7), method = method),
      f$fit$obdata, f$fit$xcoord, f$fit$ycoord, n)
    expect_identical(local$random_blocks, unname(split(seq_len(n), labels)))
    if (method == "random") expect_lte(diff(range(lengths(local$random_blocks))), 1)
    prepared <- prepare_conditional_random(f$fit, design, local)
    precision <- diag(1 / design$variance)
    for (rows in local$random_blocks) precision <- precision +
      t(Z[rows, , drop = FALSE]) %*% solve(S[rows, rows], Z[rows, , drop = FALSE])
    expect_equal(prepared$covariance, solve(precision), ignore_attr = TRUE)
    set.seed(1)
    one <- get_local_list_conditional(list(groups_random = 1, method_random = method, ordering = "none"), f$fit, f$new)
    expect_identical(one$random_blocks, list(seq_len(n)))
    expect_equal(prepare_conditional_random(f$fit, design, one)$covariance,
      prepare_conditional_random(f$fit, design)$covariance)
    set.seed(1)
    first <- get_local_list_conditional(list(groups_random = 7, method_random = method, ordering = "random"), f$fit, f$new)
    set.seed(1)
    second <- get_local_list_conditional(list(groups_random = 7, method_random = method, ordering = "random"), f$fit, f$new[1, , drop = FALSE])
    expect_identical(first$random_blocks, second$random_blocks)
  }
  set.seed(1)
  size <- get_local_list_conditional(list(size_random = 10), f$fit, f$new)
  set.seed(1)
  both <- get_local_list_conditional(list(size_random = 10, groups_random = 7), f$fit, f$new)
  expect_identical(size, both[names(size)])
  expect_identical(size$method_random, "kmeans")
  repeated <- f$fit
  repeated$obdata$cx <- repeated$obdata$cy <- 0
  expect_error(get_local_list_conditional(list(groups_random = 2), repeated, f$new), "Too few distinct")
  set.seed(1)
  random <- get_local_list_conditional(list(groups_random = 7, method_random = "random"), repeated, f$new)
  expect_length(random$random_blocks, 7)
  expect_equal(sort(unlist(random$random_blocks)), seq_len(n))
  expect_equal(lengths(conditional_random_blocks(repeated, groups = n, method = "random")), rep(1L, n))
})

test_that("effect grouping inputs reject invalid settings", {
  f <- random_conditional_fixture()
  for (name in c("size_random", "groups_random")) {
    for (value in list(0, -1, NA, Inf, 1.5, c(1, 2), "3")) {
      expect_error(get_local_list_conditional(setNames(list(value), name), f$fit, f$new), "positive integer")
    }
  }
  expect_error(get_local_list_conditional(list(groups_random = f$fit$n + 1), f$fit, f$new), "cannot exceed")
  for (method in list("none", NA, 1, c("kmeans", "random"))) {
    expect_error(get_local_list_conditional(list(method_random = method), f$fit, f$new), "method_random must")
  }
  for (setting in list(list(size_random = 10), list(groups_random = 3), list(method_random = "random"))) {
    expect_error(get_local_list_conditional(c(list(method = "all"), setting), f$fit, f$new), "only used")
    expect_error(get_local_list_conditional(c(list(approximation = "low-rank"), setting), f$fit, f$new), "only used")
    no_random <- f$fit; no_random$random <- NULL
    expect_error(get_local_list_conditional(setting, no_random, f$new), "only used")
  }
})

test_that("random grouping supports shared LM and GLM effect draws", {
  for (glm in c(FALSE, TRUE)) {
    f <- random_conditional_fixture(glm)
    for (method in c("kmeans", "random")) {
      local <- list(size = 5, groups_random = 3, method_random = method)
      set.seed(1)
      basic <- conditional(f$fit, f$new, samples = 5, local = local)
      set.seed(1)
      all <- conditional(f$fit, f$new, samples = 5, local = local, output = c("all", "random"))
      expect_identical(basic, all$newdata)
      expect_equal(dim(all$random), c(4, 5))
      expect_true(all(is.finite(all$random)))
      expect_true(all(is.finite(all$newdata)))
    }
  }
})

test_that("exact joint effect and prediction draws match Gaussian moments", {
  for (glm in c(FALSE, TRUE)) {
    f <- random_conditional_fixture(glm)
    fit <- f$fit
    design <- get_conditional_random_design(fit, f$new)
    Z <- as.matrix(design$observed)
    Zu <- as.matrix(design$newdata)
    Gamma <- diag(design$variance)
    Sigma <- as.matrix(covmatrix(fit))
    spatial <- conditional_spatial_object(fit)
    S <- as.matrix(covmatrix(spatial))
    cross <- as.matrix(covmatrix(spatial, f$new, cov_type = "obs.pred"))
    K <- t(cross) %*% solve(S)
    R <- as.matrix(covmatrix(spatial, f$new, cov_type = "pred.pred")) - K %*% cross
    T <- Gamma %*% t(Z) %*% solve(Sigma)
    H <- Gamma - T %*% Z %*% Gamma
    X <- model.matrix(fit)
    Xu <- model.matrix(~ x, f$new)
    C <- vcov(fit)
    center <- if (glm) get_latent_center(fit) else
      list(eta = model.response(model.frame(fit)), betahat = coef(fit))
    w <- center$eta
    w <- w - fit$obdata$off
    V <- matrix(0, fit$n, fit$n)
    M <- matrix(0, fit$n, NCOL(X))
    if (glm) {
      eta <- fitted(fit, type = "link")
      D <- get_D(fit$family, eta, fit$y, fit$size, as.vector(coef(fit, type = "dispersion")))
      V <- solve(solve(Sigma) - D)
      M <- V %*% solve(Sigma, X)
    }
    r <- w - X %*% center$betahat
    J <- M - X
    A <- Zu - K %*% Z
    L <- rbind(T, K + A %*% T)
    F <- rbind(T %*% J, Xu + (K + A %*% T) %*% J)
    E <- rbind(diag(NROW(Gamma)), A)
    target <- F %*% C %*% t(F) + L %*% V %*% t(L) + E %*% H %*% t(E)
    i <- NROW(Gamma) + seq_len(NROW(Xu))
    target[i, i] <- target[i, i] + R
    mu <- c(T %*% r, Xu %*% center$betahat + (K + A %*% T) %*% r + f$new$off)
    set.seed(1)
    out <- conditional(fit, f$new, samples = 6000, output = c("random", "newdata", "beta", if (glm) "latent"))
    joint <- rbind(out$random, out$newdata)
    expect_lt(max(abs(rowMeans(joint) - mu) / sqrt(diag(target))), .055)
    empirical <- cov(t(joint))
    expect_lt(max(abs(empirical - target) / sqrt(outer(diag(target), diag(target)))), .075)
    # Integrating a out recovers the original total-covariance prediction law.
    Ktotal <- t(as.matrix(covmatrix(fit, f$new, cov_type = "obs.pred"))) %*% solve(Sigma)
    Rtotal <- as.matrix(covmatrix(fit, f$new, cov_type = "pred.pred")) - Ktotal %*% Sigma %*% t(Ktotal)
    Ftotal <- Xu + Ktotal %*% J
    expect_equal(as.matrix(target[i, i]), as.matrix(Rtotal + Ktotal %*% V %*% t(Ktotal) + Ftotal %*% C %*% t(Ftotal)),
      ignore_attr = TRUE, tolerance = 1e-9)
    if (glm) {
      target_cross <- T %*% (V + J %*% C %*% t(M))
      empirical_cross <- cov(t(out$random), t(out$latent))
      sd_cross <- sqrt(outer(diag(target)[seq_len(NROW(Gamma))], diag(V + M %*% C %*% t(M))))
      expect_lt(max(abs(empirical_cross - target_cross) / sd_cross), .075)
    }
  }
})

test_that("outputs, local dispatch and zero variances behave properly", {
  for (glm in c(FALSE, TRUE)) {
    f <- random_conditional_fixture(glm)
    fit <- f$fit
    for (local in list(FALSE, list(size = 5, size_random = 10, ordering = "none"))) {
      set.seed(1)
      basic <- conditional(fit, f$new, samples = 7, local = local)
      set.seed(1)
      all <- conditional(fit, f$new, samples = 7, local = local, output = c("all", "random", if (glm) "latent"))
      expect_identical(basic, all$newdata)
      expect_identical(dim(all$random), c(4L, 7L))
      expect_identical(attr(all$random, "effect_info")$term, rep("1 | g", 4))
      expect_identical(dim(conditional(fit, samples = 1, local = local, output = "random")), c(4L, 1L))
    }
    low <- list(approximation = "low-rank", size_base = 12, size_new = 3, reorder_base = "none", reorder_new = "none")
    expect_error(conditional(fit, f$new, local = low, output = "random"), "sampled base")
    expect_true(all(is.finite(conditional(fit, f$new, samples = 2, local = low))))
    full <- list(approximation = "low-rank", method_base = "all", method_new = "all")
    expect_identical(dim(conditional(fit, f$new, samples = 1, local = full, output = "random")), c(4L, 1L))
    expect_error(conditional(fit, f$new, local = c(low, list(size_random = 10))), "only used")
    fit$coefficients$randcov[] <- 0
    expect_equal(as.numeric(conditional(fit, samples = 2, output = "random")), rep(0, 8))
    fit$random <- NULL
    fit$coefficients$randcov <- NULL
    for (local in list(FALSE, list(size = 5), low)) {
      for (output in list("random", c("beta", "random"), c("all", "random"))) {
        set.seed(1)
        rng <- .Random.seed
        expect_error(conditional(fit, samples = 2, output = output, local = local),
          'output = "random" requires a model with random effects.', fixed = TRUE)
        expect_identical(.Random.seed, rng)
      }
    }
  }
})

test_that("areal effect draws keep the full neighborhood and permit random-only output", {
  set.seed(1)
  d <- data.frame(x = seq(-1, 1, length.out = 12), g = factor(rep(1:3, 4)),
    y = c(2, NA, 4, 3, 6, 2, 1, 4, NA, 2, 3, 6))
  W <- 1 * (abs(outer(1:12, 1:12, "-")) == 1)
  for (type in c("car", "sar")) for (glm in c(FALSE, TRUE)) {
    args <- list(formula = y ~ x, data = d, W = W, random = ~ (1 | g),
      spcov_initial = spcov_initial(type, de = .4, ie = .2, range = .2, known = "given"),
      randcov_initial = randcov_initial(g = .3, known = "given"))
    if (glm) args$family <- "poisson"
    fit <- do.call(if (glm) spgautor else spautor, args)
    out <- conditional(fit, samples = 3, output = c("all", "random", if (glm) "latent"))
    expect_identical(dim(out$newdata), c(2L, 3L))
    expect_identical(rownames(out$newdata), c("2", "9"))
    expect_identical(dim(out$random), c(3L, 3L))
    expect_true(all(is.finite(out$newdata)))
    expect_error(conditional(fit, local = FALSE), "only exact")
    args$data$y[c(2, 9)] <- 3
    fit_full <- do.call(if (glm) spgautor else spautor, args)
    expect_identical(dim(conditional(fit_full, samples = 1, output = "random")), c(3L, 1L))
    fit_full$random <- NULL
    fit_full$coefficients$randcov <- NULL
    for (output in list("random", c("beta", "random"), c("all", "random"))) {
      set.seed(1)
      rng <- .Random.seed
      expect_error(conditional(fit_full, samples = 1, output = output),
        'output = "random" requires a model with random effects.', fixed = TRUE)
      expect_identical(.Random.seed, rng)
    }
  }
})

test_that("local fits, sf data and binomial transformations retain shared draws", {
  f <- random_conditional_fixture()
  d <- f$data
  d$y <- seq_len(NROW(d)) %% 5
  fit <- spglm(cbind(y, 5 - y) ~ x + offset(off), d[1:30, ], family = "binomial",
    xcoord = cx, ycoord = cy, random = ~ (1 | g), local = list(method = "kmeans", groups = 3),
    spcov_initial = spcov_initial("exponential", de = .3, ie = .2, range = .3, known = "given"),
    randcov_initial = randcov_initial(g = .3, known = "given"))
  new <- sf::st_as_sf(d[31:36, ], coords = c("cx", "cy"))
  for (local in list(FALSE, list(size = 5, size_random = 10),
      list(approximation = "low-rank", size_base = 12, size_new = 2, reorder_base = "none", reorder_new = "none"))) {
    output <- if (is.list(local) && identical(local$approximation, "low-rank")) c("newdata", "beta") else c("newdata", "beta", "random")
    set.seed(1)
    link <- conditional(fit, new, samples = 3, output = output, local = local)
    set.seed(1)
    response <- conditional(fit, new, samples = 3, output = output, local = local, type = "response")
    expect_equal(response$newdata, plogis(link$newdata))
    expect_identical(response$beta, link$beta)
    if ("random" %in% output) expect_identical(response$random, link$random)
    draws <- conditional(fit, new, samples = 3, local = local, type = "new", newdata_size = 1:6)
    expect_true(all(draws >= 0 & draws <= 1:6 & draws == floor(draws)))
  }
})

test_that("crossed slopes, nesting and new partitions keep their meaning", {
  f <- random_conditional_fixture(partition = TRUE)
  for (formula in list(~ (1 + x | g) + (1 | h), ~ (1 | g/h))) {
    fit <- f$fit
    fit$random <- formula
    names <- get_randcov_names(formula)
    fit$coefficients$randcov <- do.call(randcov_params, setNames(as.list(rep(.2, length(names))), names))
    new <- f$new
    new$part <- as.character(new$part)
    new$part[1:2] <- "unseen"
    design <- get_conditional_random_design(fit, new)
    Z <- as.matrix(design$observed)
    expect_equal(as.matrix(covmatrix(fit) - covmatrix(conditional_spatial_object(fit))),
      Z %*% diag(design$variance) %*% t(Z), ignore_attr = TRUE)
    for (local in list(FALSE, list(size = 5, size_random = 10),
        list(approximation = "low-rank", size_base = 12, size_new = 3, reorder_base = "none", reorder_new = "none"))) {
      expect_true(all(is.finite(conditional(fit, new, samples = 3, local = local))))
    }
  }
})

test_that("effect preparation is reused and covariance blocks stay local", {
  check <- function(glm) {
    f <- random_conditional_fixture(glm)
    original <- f$fit
    prepare <- prepare_conditional_random
    covfun <- if (glm) covmatrix.spglm else covmatrix.splm
    calls <- 0L
    sizes <- integer()
    bindings <- list(prepare_conditional_random = function(...) { calls <<- calls + 1L; prepare(...) })
    bindings[[if (glm) "covmatrix.spglm" else "covmatrix.splm"]] <- function(object, ...) {
      sizes <<- c(sizes, NROW(object$obdata)); covfun(object, ...)
    }
    do.call(testthat::local_mocked_bindings, c(bindings, list(.env = environment())))
      for (chunk in c(1, 20)) {
        calls <- 0L; sizes <- integer()
        conditional(f$fit, f$new, samples = 7, local = list(size = 5, size_random = 10, chunk_size = chunk))
        expect_equal(calls, 1L)
        expect_lt(max(sizes), f$fit$n)
      }
    expect_identical(f$fit, original)
    set.seed(1)
    first <- get_local_list_conditional(list(size = 5, size_random = 10, ordering = "random"), f$fit, f$new)
    set.seed(1)
    second <- get_local_list_conditional(list(size = 5, size_random = 10, ordering = "random"), f$fit, f$new[1, , drop = FALSE])
    expect_identical(first$random_blocks, second$random_blocks)
    low <- list(approximation = "low-rank", size_base = 12, size_new = 2,
      reorder_base = "none", reorder_new = "none", chunk_size = 2)
    set.seed(1)
    serial <- conditional(f$fit, f$new, samples = 7, local = low)
    set.seed(1)
    parallel <- conditional(f$fit, f$new, samples = 7, local = c(low, list(parallel = TRUE, ncores = 2)))
    expect_equal(serial, parallel)
  }
  check(FALSE)
  check(TRUE)
})

test_that("theta draws use their own effects and retain randcov parameter output", {
  f <- random_conditional_fixture()
  fit <- f$fit
  # A deterministic theta sequence isolates rebuilding from estimation.
  fit$is_known$randcov[] <- FALSE
  fit$vcov$cov <- matrix(.01, 1, 1, dimnames = list("1 | g", "1 | g"))
  prepare <- prepare_conditional_random
  variances <- list()
  counter <- 0L
  local_mocked_bindings(simulate_theta_draw = function(...) {
    counter <<- counter + 1L
    list(spcov_params_b = coef(fit, type = "spcov"),
      randcov_params_b = randcov_params("1 | g" = if (counter %% 2L) 0 else .7))
  }, prepare_conditional_random = function(object, design, ...) {
    variances[[length(variances) + 1L]] <<- design$variance
    prepare(object, design, ...)
  })
  set.seed(1)
  out <- conditional(fit, f$new, samples = 4, simulate_covparams = TRUE,
    output = c("random", "newdata", "randcov", "cov", "beta"))
  expect_equal(as.numeric(out$randcov), c(0, .7, 0, .7))
  expect_equal(as.numeric(out$random[, c(1, 3)]), rep(0, 8))
  expect_equal(variances, lapply(c(0, .7, 0, .7), rep, 4))
  expect_true(all(is.finite(out$newdata)))
  expect_equal(unname(out$cov["1 | g", ]), as.numeric(out$randcov))
  expect_identical(dim(conditional(fit, output = "random", samples = 1, simulate_covparams = TRUE)), c(4L, 1L))
})

test_that("fixed theta composition agrees with analytical LM prediction moments", {
  f <- random_conditional_fixture()
  fit <- splm(y ~ x + offset(off), f$data[1:12, ], xcoord = cx, ycoord = cy,
    random = ~ (1 | g),
    spcov_initial = spcov_initial("exponential", de = .3, ie = .2, range = .3, known = "given"),
    randcov_initial = randcov_initial(g = .3, known = "given"))
  fit$is_known$randcov[] <- FALSE
  fit$vcov$cov <- matrix(.01, 1, 1, dimnames = list("1 | g", "1 | g"))
  local_mocked_bindings(simulate_theta_draw = function(...) list(
    spcov_params_b = coef(fit, type = "spcov"), randcov_params_b = coef(fit, type = "randcov")))
  set.seed(1)
  out <- conditional(fit, f$new, samples = 1500, simulate_covparams = TRUE,
    output = c("newdata", "beta", "random"))
  target <- predict(fit, f$new, se.fit = TRUE, local = FALSE)
  expect_lt(max(abs(rowMeans(out$newdata) - target$fit) / target$se.fit), .1)
  expect_lt(max(abs(apply(out$newdata, 1, sd) / target$se.fit - 1)), .075)
  expect_lt(max(abs(cov(t(out$beta)) - vcov(fit)) / sqrt(outer(diag(vcov(fit)), diag(vcov(fit))))), .1)
})
