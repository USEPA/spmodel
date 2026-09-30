skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

joint_fixture <- function(local = FALSE, random = FALSE, family = "poisson", partition = FALSE,
                          coordinates = "data", intercept = FALSE, missing = FALSE, near = FALSE) {
  data <- data.frame(cx = seq(0, 1, length.out = 18),
    cy = rep(c(0, 0.3, 0.1), 6), x = rep(c(-1, 0, 1), 6),
    off = seq(-0.3, 0.2, length.out = 18), group = factor(rep(1:3, 6)),
    y = c(2, 1, 4, 3, 6, 2, 1, 4, 5, 2, 3, 6, 4, 2, 5, 3, 1, 7))
  if (family == "binomial") data$y <- data$y %% 5
  if (family == "beta") data$y <- data$y / 10
  newdata <- data[c(2, 6, 11, 15), ]
  newdata$cx <- newdata$cx + 0.025
  data$part <- factor(rep(c(1, 2), 9), levels = 1:3)
  newdata$part <- data$part[c(2, 6, 11, 15)]
  if (missing) data$y[c(2, 6)] <- NA
  if (coordinates == "one") {
    data$cy <- 0
    newdata$cy <- 0
  }
  if (coordinates == "sf") {
    data <- sf::st_as_sf(data, coords = c("cx", "cy"))
    newdata <- sf::st_as_sf(newdata, coords = c("cx", "cy"))
  }
  formula <- if (family == "binomial") cbind(y, 5 - y) ~ x + offset(off) else y ~ x + offset(off)
  if (intercept) formula <- y ~ 1 + offset(off)
  args <- list(formula = formula, data = data, local = local,
    spcov_initial = spcov_initial("exponential", de = 0.4, ie = if (near) 1e-8 else 0.2,
      range = if (near) 1e4 else 0.3, known = "given"),
    partition_factor = if (partition) ~ part else NULL,
    random = if (random) ~ (1 | group) + (0 + x | group) else NULL,
    randcov_initial = if (random) randcov_initial(group = 0.2, "0 + x | group" = 0.1, known = "given") else NULL,
    dispersion_initial = dispersion_initial(as.character(family),
      dispersion = if (family %in% c("poisson", "binomial")) 1 else 3, known = "given"))
  if (coordinates != "sf") args$xcoord <- "cx"
  if (coordinates == "data") args$ycoord <- "cy"
  fit <- do.call(spglm, args)
  list(fit = fit, newdata = newdata)
}

joint_reference <- function(fit, newdata, base = seq_len(fit$n), blocks = list(seq_len(NROW(newdata))), neighbors = NULL) {
  X <- model.matrix(fit)
  Xnew <- model.matrix(delete.response(terms(fit)), newdata, contrasts.arg = fit$contrasts)
  Sigma <- as.matrix(covmatrix(fit))
  P <- matrix(0, fit$n, fit$n)
  index <- fit$local_index
  if (is.null(index) && !is.null(fit$partition_factor)) index <- model.frame(fit$partition_factor, fit$obdata)[[1L]]
  groups <- if (is.null(index)) list(seq_len(fit$n)) else split(seq_len(fit$n), index, drop = TRUE)
  for (rows in groups) P[rows, rows] <- solve(Sigma[rows, rows])
  information <- t(X) %*% P %*% X
  if (length(groups) > 1L) information <- information + diag(fit$diagtol, NCOL(X))
  G <- solve(information)
  Vu <- vcov(fit, var_correct = FALSE)
  B <- Vu %*% t(X) %*% P
  L <- solve(P - P %*% X %*% G %*% t(X) %*% P + diag(exp(fitted(fit, type = "link"))))
  C <- as.matrix(covmatrix(fit, newdata))
  K <- as.matrix(covmatrix(fit, newdata, cov_type = "pred.pred"))
  m <- NROW(newdata)
  A <- matrix(0, m, fit$n)
  R <- matrix(0, m, m)
  if (is.null(neighbors)) {
    A[, base] <- C[, base, drop = FALSE] %*% solve(Sigma[base, base, drop = FALSE])
    for (rows in blocks) R[rows, rows] <- K[rows, rows, drop = FALSE] - A[rows, , drop = FALSE] %*% t(C[rows, , drop = FALSE])
  } else {
    joint_cov <- rbind(cbind(Sigma, t(C)), cbind(C, K))
    coords <- rbind(fit$obdata[c("cx", "cy")], newdata[c("cx", "cy")])
    T <- matrix(0, m, m)
    innovations <- numeric(m)
    for (i in seq_len(m)) {
      pool <- seq_len(fit$n + i - 1L)
      d <- rowSums((as.matrix(coords[pool, ]) - matrix(as.numeric(coords[fit$n + i, ]), length(pool), 2, byrow = TRUE))^2)
      pool <- pool[order(d)[seq_len(min(neighbors, length(pool)))]]
      weight <- solve(joint_cov[pool, pool], joint_cov[pool, fit$n + i])
      observed <- pool <= fit$n
      A[i, pool[observed]] <- weight[observed]
      T[i, pool[!observed] - fit$n] <- weight[!observed]
      innovations[i] <- K[i, i] - sum(weight * joint_cov[pool, fit$n + i])
    }
    propagation <- solve(diag(m) - T)
    A <- propagation %*% A
    R <- propagation %*% diag(innovations) %*% t(propagation)
  }
  H <- Xnew - A %*% X
  W <- H %*% B + A
  beta_cov <- Vu + B %*% L %*% t(B)
  pred_cov <- R + H %*% Vu %*% t(H) + W %*% L %*% t(W)
  cross <- Vu %*% t(H) + B %*% L %*% t(W)
  w <- fitted(fit, type = "link") - fit$obdata$off
  mu <- Xnew %*% coef(fit) + A %*% (w - X %*% coef(fit)) + newdata$off
  list(P = P, G = G, B = B, L = L, Vu = Vu, H = H,
    mean = c(coef(fit), mu), covariance = rbind(cbind(beta_cov, cross), cbind(t(cross), pred_cov)))
}

expect_joint_moments <- function(draws, target) {
  s <- NCOL(draws)
  S <- target$covariance
  mcse <- sqrt((S^2 + outer(diag(S), diag(S))) / (s - 1))
  expect_lt(max(abs(cov(t(draws)) - S) / mcse), 6)
  expect_lt(max(abs(rowMeans(draws) - target$mean) / sqrt(diag(S) / s)), 6)
}

test_that("exact joint factors and full sampled moments match independent matrices", {
  fixture <- joint_fixture()
  fit <- fixture$fit
  target <- joint_reference(fit, fixture$newdata)
  joint <- get_conditional_glm_joint(fit)
  expect_equal(unname(joint$wts_beta), unname(target$B), tolerance = 1e-10)
  expect_equal(chol2inv(t(joint$cov_lowchol_mH)), target$L, tolerance = 1e-10)
  expect_equal(unname(vcov(fit)), unname(target$covariance[1:2, 1:2]), tolerance = 1e-4)
  prediction <- predict(fit, fixture$newdata, se.fit = TRUE, type = "link", local = FALSE)
  expect_equal(as.numeric(prediction$fit), unname(target$mean[-(1:2)]), tolerance = 1e-8)
  expect_equal(as.numeric(prediction$se.fit)^2, unname(diag(target$covariance)[-(1:2)]), tolerance = 1e-8)
  excess <- target$H %*% target$B %*% target$L %*% t(target$B) %*% t(target$H)
  expect_gt(max(abs(excess)), 0.001)
  set.seed(832)
  draws <- conditional(fit, fixture$newdata, output = "all", samples = 30000, local = FALSE)
  expect_joint_moments(rbind(draws$beta, draws$newdata), target)
  expect_equal(draws$object[, 1], unname(fitted(fit, type = "link")))
  expect_identical(rownames(draws$beta), names(coef(fit)))
})

test_that("local fitting and simulation preserve the constructed joint moments", {
  simulation <- list(FALSE,
    list(method_base = "all", method_new = "all"),
    list(method_base = "base", size_base = 9, reorder_base = "none", method_new = "base", size_new = 2, reorder_new = "none"),
    list(approximation = "vecchia", method = "all", ordering = "none"),
    list(approximation = "vecchia", method = "distance", size = 5, ordering = "none"))
  for (adjustment in c("exact", "none", "theoretical", "pooled")) {
    local_fit <- if (adjustment == "exact") FALSE else list(index = rep(c(3, 1, 2, 1, 2, 3), 3), var_adjust = adjustment)
    fixture <- joint_fixture(local_fit)
    fit <- fixture$fit
    for (i in seq_along(simulation)) {
      target <- joint_reference(fit, fixture$newdata,
        base = if (i == 3) 1:9 else 1:18,
        blocks = if (i == 3) list(1:2, 3:4) else list(1:4),
        neighbors = if (i == 5) 5 else NULL)
      joint <- get_conditional_glm_joint(fit)
      expect_equal(chol2inv(t(joint$cov_lowchol_mH)), target$L, tolerance = 1e-9)
      expect_equal(unname(joint$wts_beta), unname(target$B), tolerance = 1e-9)
      expect_equal(unname(vcov(fit)), unname(target$covariance[1:2, 1:2]), tolerance = 1e-4)
      set.seed(500 + i)
      draws <- conditional(fit, fixture$newdata, output = c("beta", "newdata"), samples = 30000, local = simulation[[i]])
      expect_joint_moments(rbind(draws$beta, draws$newdata), target)
    }
  }
})

test_that("coupled draws preserve family, offset, dimension and snapshots", {
  for (family in c("poisson", "nbinomial", "binomial", "Gamma", "inverse.gaussian", "beta")) {
    fixture <- joint_fixture(family = family)
    for (local in list(FALSE, list(approximation = "vecchia", size = 5, ordering = "none"))) {
      fit <- fixture$fit
      new <- fixture$newdata[1, , drop = FALSE]
      set.seed(82)
      link <- conditional(fit, new, samples = 1, local = local, output = "all")
      set.seed(82)
      response <- conditional(fit, new, samples = 1, local = local, type = "response", newdata_size = 5)
      expect_equal(dim(link$newdata), c(1L, 1L))
      expect_equal(dim(response), c(1L, 1L))
      expect_equal(as.numeric(response), as.numeric(invlink(link$newdata, family, size = 5)))
      expect_equal(link$object[, 1], unname(fitted(fit, type = "link")))
      set.seed(82)
      expect_identical(conditional(fit, new, samples = 1, local = local, output = "all"), link)
      set.seed(82)
      observation <- conditional(fit, new, samples = 1, local = local, type = "new", newdata_size = 5)
      expect_true(all(is.finite(observation)))
      if (family %in% c("poisson", "binomial", "nbinomial")) expect_equal(observation, round(observation))
      if (family == "binomial") expect_true(all(observation >= 0 & observation <= 5))
      if (family == "beta") expect_true(all(observation > 0 & observation < 1))
    }
  }
})

test_that("random components, partitions and parallel blocks retain shared uncertainty", {
  fixture <- joint_fixture(random = TRUE, partition = TRUE)
  fit <- fixture$fit
  original <- fit
  for (simulation in list(FALSE,
      list(method_base = "base", size_base = 9, reorder_base = "none",
        method_new = "base", size_new = 2, reorder_new = "none", parallel = TRUE, ncores = 2),
      list(approximation = "vecchia", method = "distance", size = 5, ordering = "none"))) {
    lowrank <- is.list(simulation) && isTRUE(simulation$parallel)
    vecchia <- is.list(simulation) && identical(simulation$approximation, "vecchia")
    target <- joint_reference(fit, fixture$newdata, base = if (lowrank) 1:9 else 1:18,
      blocks = if (lowrank) list(1:2, 3:4) else list(1:4), neighbors = if (vecchia) 5 else NULL)
    set.seed(489)
    draws <- conditional(fit, fixture$newdata, output = "all", samples = 30000, local = simulation)
    expect_joint_moments(rbind(draws$beta, draws$newdata), target)
    expect_identical(fit, original)
  }
})

test_that("coordinate forms, missing responses and intercept-only fits remain supported", {
  for (coordinates in c("data", "sf", "one")) {
    fixture <- joint_fixture(coordinates = coordinates, intercept = TRUE, missing = TRUE)
    fit <- fixture$fit
    for (local in list(FALSE, list(approximation = "vecchia", method = "all", ordering = "none"))) {
      set.seed(141)
      explicit <- conditional(fit, fit$newdata, output = "all", samples = 1, local = local)
      set.seed(141)
      implicit <- conditional(fit, output = "all", samples = 1, local = local)
      expect_identical(explicit, implicit)
      expect_equal(dim(explicit$beta), c(1L, 1L))
      expect_equal(dim(explicit$newdata), c(2L, 1L))
      expect_equal(explicit$object[, 1], unname(fitted(fit, type = "link")))
      expect_true(all(is.finite(explicit$newdata)))
    }
  }
})

test_that("near-singular spatial covariance and invalid latent precision are explicit", {
  fixture <- joint_fixture(near = TRUE)
  expect_true(all(is.finite(conditional(fixture$fit, fixture$newdata, samples = 20))))
  fixture <- joint_fixture()
  testthat::local_mocked_bindings(get_D = function(...) diag(1e6, 18))
  expect_error(get_conditional_glm_joint(fixture$fit), "latent precision is not positive definite")
})

test_that("factor contrasts and parallel draws are reproducible", {
  fixture <- joint_fixture()
  data <- fixture$fit$obdata
  data$x <- seq(-1, 1, length.out = NROW(data))
  fit <- spglm(y ~ group + x + offset(off), family = "poisson", data = data,
    xcoord = "cx", ycoord = "cy", contrasts = list(group = "contr.sum"),
    spcov_initial = spcov_initial("exponential", de = 0.4, ie = 0.2, range = 0.3, known = "given"))
  newdata <- fixture$newdata[c(3, 2, 1, 4), ]
  set.seed(382)
  draws <- conditional(fit, newdata, output = "all", samples = 30000, local = FALSE)
  target <- joint_reference(fit, newdata)
  expect_joint_moments(rbind(draws$beta, draws$newdata), target)
  expect_identical(rownames(draws$beta), names(coef(fit)))
  local <- list(method_base = "base", size_base = 9, reorder_base = "random",
    method_new = "base", size_new = 2, reorder_new = "random", kmeans_new = FALSE,
    parallel = TRUE, ncores = 2)
  set.seed(834)
  first <- conditional(fit, newdata, output = "all", samples = 10, local = local)
  set.seed(834)
  expect_identical(conditional(fit, newdata, output = "all", samples = 10, local = local), first)
})
