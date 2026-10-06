skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

# Extra arguments are only defined when extras are enabled.

areal_conditional_fixture <- function(family = NULL, spcov_type = "car", row_st = TRUE,
                                      missing = c(2, 7, 15), formula = y ~ x + offset(off),
                                      random = FALSE, partition = FALSE, polygons = FALSE) {
  n <- 18L
  data <- data.frame(x = seq(-1, 1, length.out = n), off = seq(-0.2, 0.3, length.out = n),
    group = factor(rep(1:3, 6)), part = factor(rep(1:2, each = 9)),
    y = c(2, 1, 4, 3, 6, 2, 1, 4, 5, 2, 3, 6, 4, 2, 5, 3, 1, 7))
  if (identical(family, "binomial")) {
    data$y <- data$y %% 4
    data$failures <- 5 - data$y
    formula <- cbind(y, failures) ~ x + offset(off)
  }
  if (identical(family, "beta")) data$y <- data$y / 10
  data$y[missing] <- NA
  if (identical(family, "binomial") && length(missing) > 1L) {
    data$y[missing[2L]] <- 2
    data$failures[missing[2L]] <- NA
  }
  W <- 1 * (abs(outer(seq_len(n), seq_len(n), "-")) == 1)
  if (polygons) {
    geometry <- sf::st_make_grid(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 6, ymax = 3)), n = c(6, 3))
    data <- sf::st_sf(data, geometry = geometry)
  }
  args <- list(formula = formula, data = data, W = W, row_st = row_st,
    spcov_initial = spcov_initial(spcov_type, de = 0.4, ie = 0.2, range = 0.2, known = "given"))
  if (!row_st && spcov_type == "car") args$M <- rep(1, n)
  if ("group" %in% all.vars(formula)) args$contrasts <- list(group = "contr.sum")
  if (random) {
    args$random <- ~ (1 | group) + (0 + x | group)
    args$randcov_initial <- randcov_initial(group = 0.2, "0 + x | group" = 0.1, known = "given")
  }
  if (partition) args$partition_factor <- ~ part
  if (is.null(family)) {
    args$ddf <- "asymptotic"
    do.call(spautor, args)
  } else {
    args$dispersion_initial <- dispersion_initial(as.character(family),
      dispersion = if (family %in% c("poisson", "binomial")) 1 else 3, known = "given")
    do.call(spgautor, args)
  }
}

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

joint_reference <- function(fit, newdata, base = seq_len(fit$n), blocks = list(seq_len(NROW(newdata))), neighbors = NULL, base_latent = FALSE) {
  if (base_latent) return(lowrank_joint_reference(fit, newdata, base, blocks))
  if (!is.null(neighbors)) return(vecchia_joint_reference(fit, newdata, size = neighbors))
  X <- model.matrix(fit)
  Xnew <- model.matrix(delete.response(terms(fit)), newdata, contrasts.arg = fit$contrasts)
  Sigma <- as.matrix(covmatrix(fit))
  SigInv <- solve(Sigma)
  eta <- fitted(fit, type = "link")
  curvature <- switch(fit$family, poisson = exp(eta),
    binomial = rowSums(model.response(model.frame(fit))) * plogis(eta) * (1 - plogis(eta)))
  H <- solve(SigInv + diag(curvature, NROW(X)))
  M <- H %*% SigInv %*% X
  C <- vcov(fit)
  cross <- as.matrix(covmatrix(fit, newdata))
  A <- cross %*% SigInv
  R <- as.matrix(covmatrix(fit, newdata, cov_type = "pred.pred")) - A %*% t(cross)
  F <- Xnew - A %*% X + A %*% M
  pred_cov <- R + A %*% H %*% t(A) + F %*% C %*% t(F)
  beta_cross <- C %*% t(F)
  offset <- model.offset(model.frame(fit))
  w <- eta - if (is.null(offset)) 0 else offset
  offset_new <- model.offset(model.frame(delete.response(terms(fit)), newdata))
  mu <- as.numeric(Xnew %*% coef(fit) + A %*% (w - X %*% coef(fit))) +
    if (is.null(offset_new)) 0 else offset_new
  list(H = H, M = M, C = C, F = F,
    mean = c(coef(fit), mu), covariance = rbind(cbind(C, beta_cross),
      cbind(t(beta_cross), pred_cov)))
}

lowrank_joint_reference <- function(fit, newdata, base, blocks = list(seq_len(NROW(newdata)))) {
  X <- model.matrix(fit)[base, , drop = FALSE]
  Xnew <- model.matrix(delete.response(terms(fit)), newdata, contrasts.arg = fit$contrasts)
  Sigma <- as.matrix(covmatrix(fit))[base, base, drop = FALSE]
  cross <- as.matrix(covmatrix(fit, newdata))[, base, drop = FALSE]
  K <- as.matrix(covmatrix(fit, newdata, cov_type = "pred.pred"))
  eta <- fitted(fit, type = "link")
  offset <- model.offset(model.frame(fit))
  w <- (eta - if (is.null(offset)) 0 else offset)[base]
  curvature <- switch(fit$family, poisson = exp(eta),
    binomial = rowSums(model.response(model.frame(fit))) * plogis(eta) * (1 - plogis(eta)))
  P <- solve(Sigma)
  H <- solve(P + diag(curvature[base], length(base)))
  M <- H %*% P %*% X
  V <- vcov(fit)
  A <- cross %*% P
  F <- Xnew - A %*% X + A %*% M
  R <- matrix(0, NROW(newdata), NROW(newdata))
  for (rows in blocks) R[rows, rows] <- K[rows, rows, drop = FALSE] -
    A[rows, , drop = FALSE] %*% t(cross[rows, , drop = FALSE])
  pred_cov <- R + A %*% H %*% t(A) + F %*% V %*% t(F)
  beta_cross <- V %*% t(F)
  offset_new <- model.offset(model.frame(delete.response(terms(fit)), newdata))
  mu <- as.numeric(Xnew %*% coef(fit) + A %*% (w - X %*% coef(fit))) +
    if (is.null(offset_new)) 0 else offset_new
  list(mean = c(coef(fit), mu), covariance = rbind(cbind(V, beta_cross),
    cbind(t(beta_cross), pred_cov)), H = H, M = M)
}

expect_joint_moments <- function(draws, target) {
  s <- NCOL(draws)
  S <- target$covariance
  mcse <- sqrt((S^2 + outer(diag(S), diag(S))) / (s - 1))
  expect_lt(max(abs(cov(t(draws)) - S) / mcse), 6)
  expect_lt(max(abs(rowMeans(draws) - target$mean) / sqrt(diag(S) / s)), 6)
}

# Independent small dense reference for tests, copied into extras file.
vecchia_joint_reference <- function(fit, newdata, size = Inf, ord = seq_len(NROW(newdata)), method = "distance",
                                    order_o = seq_len(fit$n)) {
  n <- fit$n; m <- NROW(newdata); p <- length(coef(fit))
  X <- model.matrix(fit)
  Xnew <- model.matrix(delete.response(terms(fit)), newdata, contrasts.arg = fit$contrasts)[, colnames(X), drop = FALSE]
  Sigma <- covmatrix(fit); cross <- covmatrix(fit, newdata)
  full <- rbind(cbind(Sigma, t(cross)), cbind(cross, covmatrix(fit, newdata, cov_type = "pred.pred")))
  coords <- as.matrix(rbind(fit$obdata[, c(fit$xcoord, fit$ycoord)], newdata[, c(fit$xcoord, fit$ycoord)]))
  if (fit$anisotropy) {
    params <- coef(fit, type = "spcov"); angle <- params[["rotate"]]
    coords <- coords %*% matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2)
    coords[, 2] <- coords[, 2] / params[["scale"]]
  }
  partition <- if (!is.null(fit$partition_factor)) {
    columns <- all.vars(fit$partition_factor)
    as.matrix(partition_matrix(fit$partition_factor,rbind(fit$obdata[,columns,drop=FALSE],newdata[,columns,drop=FALSE])))
  } else NULL
  select <- function(i, pool, budget) {
    if (!is.null(partition)) pool <- pool[partition[i,pool]==1]
    if (!length(pool) || budget == 0) return(integer())
    distance <- sqrt(rowSums(sweep(coords[pool, , drop = FALSE], 2, coords[i, ], "-")^2))
    score <- if (method == "covariance") -abs(as.numeric(full[i,pool])) else distance
    pool[order(score, pool)[seq_len(min(budget, length(pool)))]]
  }
  eta <- fitted(fit, type = "link")
  w <- as.numeric(w_offset_free(eta, model.offset(model.frame(fit))))
  D <- diag(get_D(fit$family, eta, fit$y, fit$size, as.vector(coef(fit, type = "dispersion"))))
  Mo <- matrix(0, n, p); Lo <- matrix(0, n, n)
  for (k in seq_len(n)) {
    i <- order_o[k]; N <- select(i, order_o[seq_len(k - 1L)], size)
    K <- c(i, select(i, setdiff(seq_len(n), i), size - 1))
    A <- unique(c(i, N, K)); nn <- match(N, A)
    inv <- solve(Sigma[A, A, drop = FALSE])
    V <- solve(inv - diag(D[A], length(A))); M <- V %*% inv %*% X[A, , drop = FALSE]
    g <- if (length(N)) as.numeric(V[1, nn, drop = FALSE] %*% solve(V[nn, nn, drop = FALSE])) else numeric()
    Mo[i, ] <- M[1, ] - colSums(g * M[nn, , drop = FALSE]) + colSums(g * Mo[N, , drop = FALSE])
    Lo[i, ] <- colSums(g * Lo[N, , drop = FALSE])
    Lo[i, i] <- sqrt(V[1, 1] - sum(g * V[nn, 1]))
  }
  A <- rbind(diag(n), matrix(0, m, n)); L <- matrix(0, n + m, m)
  for (k in seq_len(m)) {
    i <- n + ord[k]; N <- select(i, c(seq_len(n), n + ord[seq_len(k - 1L)]), size)
    g <- if (length(N)) solve(full[N, N, drop = FALSE], full[N, i]) else numeric()
    A[i, ] <- as.numeric(crossprod(g, A[N, , drop = FALSE]))
    L[i, ] <- as.numeric(crossprod(g, L[N, , drop = FALSE]))
    L[i, ord[k]] <- sqrt(full[i, i] - sum(g * full[N, i]))
  }
  A <- A[n + seq_len(m), , drop = FALSE]; L <- L[n + seq_len(m), , drop = FALSE]
  F <- Xnew + A %*% (Mo - X); C <- vcov(fit)
  pred_cov <- F %*% C %*% t(F) + A %*% tcrossprod(Lo) %*% t(A) + tcrossprod(L)
  mu <- as.numeric(Xnew %*% coef(fit) + A %*% (w - X %*% coef(fit)))
  offset <- model.offset(model.frame(delete.response(terms(fit)), newdata))
  if (!is.null(offset)) mu <- mu + offset
  list(mean = c(coef(fit), mu), covariance = rbind(cbind(C, C %*% t(F)), cbind(F %*% C, pred_cov)),
    observed_mean = w, observed_covariance = Mo %*% C %*% t(Mo) + tcrossprod(Lo))
}

# General conditional simulation checks

set.seed(1)

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

exdata_pois <- exdata
exdata_pois$count <- round(abs(exdata_pois$y) * 3)

test_that("areal conditional simulation requires missing-response locations", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  spmod_autor <- spautor(y ~ x, exdata_poly, spcov_type = "car")
  expect_error(conditional(spmod_autor), "No missing data to simulate")
})

test_that("conditional() validates output/type and the newdata/object$newdata requirement", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_error(conditional(spmod, newdata = newexdata, output = "bogus"), "output must be")

  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_error(conditional(spmod_g, newdata = newexdata, output = "bogus"), "output must be")
  expect_error(conditional(spmod_g, newdata = newexdata, type = "bogus"), "should be one of")

  # no newdata argument and no missing observations to fall back on
  expect_error(conditional(spmod), "No missing data to predict")
})

test_that("output argument returns the documented shapes for splm() and spglm()", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  samples <- 20

  cond_newdata <- conditional(spmod, newdata = newexdata, output = "newdata", samples = samples)
  expect_equal(dim(cond_newdata), c(NROW(newexdata), samples))

  cond_beta <- conditional(spmod, newdata = newexdata, output = "beta", samples = samples)
  expect_equal(dim(cond_beta), c(length(coef(spmod)), samples))

  cond_object <- conditional(spmod, newdata = newexdata, output = "object", samples = samples)
  expect_equal(dim(cond_object), c(spmod$n, samples))
  # output = "object" just replicates the observed response, so every column
  # (simulation draw) is identical and equal to the observed y
  expect_true(all(apply(cond_object, 1, function(row) length(unique(row)) == 1)))
  expect_equal(cond_object[, 1], unname(model.response(model.frame(spmod))))

  cond_multi <- conditional(spmod, newdata = newexdata, output = c("newdata", "beta"), samples = samples)
  expect_named(cond_multi, c("newdata", "beta"))

  cond_all <- conditional(spmod, newdata = newexdata, output = "all", samples = samples)
  expect_named(cond_all, c("newdata", "beta", "object"))

  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  cond_object_g <- conditional(spmod_g, newdata = newexdata, output = "object", samples = samples)
  expect_equal(dim(cond_object_g), c(spmod_g$n, samples))
  # for spglm(), "object" replicates the fitted link-scale latent process w,
  # not the observed response y
  expect_equal(cond_object_g[, 1], unname(fitted(spmod_g, type = "link")))
})

test_that("type argument controls the returned scale for spglm()", {
  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  samples <- 30

  # the type-scale conversion is the last step applied to an otherwise
  # identical sequence of random draws, so matching seeds makes the
  # link/response comparison below an exact check, not just approximate
  set.seed(1)
  cond_link <- conditional(spmod_g, newdata = newexdata, type = "link", samples = samples)
  set.seed(1)
  cond_response <- conditional(spmod_g, newdata = newexdata, type = "response", samples = samples)
  cond_new <- conditional(spmod_g, newdata = newexdata, type = "new", samples = samples)

  expect_true(all(cond_response >= 0))
  expect_true(all(cond_new >= 0))
  # a "new" draw for a Poisson response is a simulated count
  expect_equal(cond_new, round(cond_new))
  # the response-scale mean is exp() of the link-scale draw here (log link,
  # no dispersion adjustment for a Poisson response)
  expect_equal(cond_response, exp(cond_link))
})

test_that("conditional.spglm() newdata_size controls the size of simulated binomial draws", {
  exdata_bin <- exdata
  exdata_bin$bern <- rbinom(NROW(exdata_bin), size = 1, prob = 0.5)
  spmod_bin <- spglm(bern ~ x, family = binomial, data = exdata_bin, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential")

  cond_default <- conditional(spmod_bin, newdata = newexdata, type = "new", samples = 30)
  expect_true(all(cond_default %in% c(0, 1)))

  cond_size5 <- conditional(spmod_bin, newdata = newexdata, type = "new", samples = 30, newdata_size = rep(5, NROW(newexdata)))
  expect_true(all(cond_size5 >= 0 & cond_size5 <= 5))
  expect_true(any(cond_size5 > 1)) # exercising sizes above the default of 1
})

test_that("conditional() respects offset() identically to an equivalent shifted-response model", {
  exdata_off <- exdata
  exdata_off$offset <- 2
  newexdata_off <- newexdata
  newexdata_off$offset <- 2

  spmod_off <- splm(y ~ x + offset(offset), exdata_off, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  spmod_no_off <- splm(I(y - 2) ~ x, exdata_off, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  set.seed(1)
  cond_off <- conditional(spmod_off, newdata = newexdata_off, samples = 50)
  set.seed(1)
  cond_no_off <- conditional(spmod_no_off, newdata = newexdata_off, samples = 50)
  # matched seeds + an additive offset should reproduce the same draws shifted
  # by exactly the offset (2), since the offset is added back on deterministically
  expect_equal(cond_off, cond_no_off + 2)
})

test_that("conditional() works with random effects, anisotropy, a partition factor, and a pure-nugget covariance", {
  samples <- 20

  spmod_rand <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~group)
  expect_error(cond_rand <- conditional(spmod_rand, newdata = newexdata, samples = samples), NA)
  expect_true(all(is.finite(cond_rand)))

  spmod_anis <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, anisotropy = TRUE)
  expect_error(cond_anis <- conditional(spmod_anis, newdata = newexdata, samples = samples), NA)
  expect_true(all(is.finite(cond_anis)))

  spmod_pf <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, partition_factor = ~group)
  expect_error(cond_pf <- conditional(spmod_pf, newdata = newexdata, samples = samples), NA)
  expect_true(all(is.finite(cond_pf)))

  # spcov_type = "none": a pure-nugget model, exercised because
  # get_conditional_new_from_base_adjust() takes a diagonal-matrix shortcut
  # (element-wise sqrt() instead of chol()) in this case
  spmod_nugget <- splm(y ~ x, exdata, spcov_type = "none", xcoord = xcoord, ycoord = ycoord)
  expect_error(cond_nugget <- conditional(spmod_nugget, newdata = newexdata, samples = samples), NA)
  expect_true(all(is.finite(cond_nugget)))
})

test_that("explicit low-rank settings run the block-processing path for splm() and spglm()", {
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  # size_base/size_new well below n = 100 / n_pred = 10 forces the block path
  # even on this small fixture (see get_local_list_conditional())
  local_small <- list(approximation = "low-rank", size_base = 30, size_new = 3)

  cond_local <- conditional(spmod, newdata = newexdata, samples = 20, local = local_small)
  expect_equal(dim(cond_local), c(NROW(newexdata), 20))
  expect_false(anyNA(cond_local))

  # base-only subsetting (method_new stays "all" since size_new >= n_pred):
  # a regression test for a fixed bug where local$index was previously only
  # set inside the method_new != "all" branch, leaving it NULL whenever a
  # large observed sample needed subsetting but a small newdata did not
  cond_local_base_only <- conditional(spmod, newdata = newexdata, samples = 20, local = list(approximation = "low-rank", size_base = 30, size_new = 500))
  expect_equal(dim(cond_local_base_only), c(NROW(newexdata), 20))
  expect_false(anyNA(cond_local_base_only))

  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  cond_local_g <- conditional(spmod_g, newdata = newexdata, samples = 20, local = local_small)
  expect_equal(dim(cond_local_g), c(NROW(newexdata), 20))
  expect_false(anyNA(cond_local_g))
})

test_that("local kmeans partitioning correctly forwards extra list elements like parallel/ncores to kmeans() (regression test)", {
  # regression test for a bug in get_local_list.R where the extra local list
  # elements left over after removing the recognized names (e.g. "parallel",
  # "ncores") were passed to kmeans() as bare names (a character vector)
  # instead of their values (local[names]); do.call() then supplied those
  # name strings as unnamed positional arguments, landing in kmeans()'s
  # nstart/algorithm slots and crashing with a match.arg() error as soon as
  # kmeans-based block partitioning ran alongside parallel = TRUE
  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_error(
    conditional(spmod_g, newdata = newexdata, samples = 20, local = list(approximation = "low-rank", size_base = 30, size_new = 3, parallel = TRUE, ncores = 2)),
    NA
  )
})

test_that("conditional() SD is close to predict() se.fit under ordinary (non-extreme) conditions", {
  # a broad sanity check, not a tight regression test (see the dedicated
  # regression tests below for that): conditional() draws should have
  # roughly the same spread as predict()'s analytic standard error
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  preds <- predict(spmod, newdata = newexdata, se.fit = TRUE)
  set.seed(1)
  cond <- conditional(spmod, newdata = newexdata, samples = 3000)
  ratio <- apply(cond, 1, sd) / preds$se.fit
  expect_true(all(ratio > 0.8 & ratio < 1.25))
  # simulated means should track predict()'s fitted values (no bias); an
  # absolute threshold is used since the fitted values straddle zero, which
  # would make a relative tolerance behave inconsistently across rows
  expect_true(all(abs(apply(cond, 1, mean) - preds$fit) < 0.1))

  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  preds_g <- predict(spmod_g, newdata = newexdata, type = "link", se.fit = TRUE)
  set.seed(1)
  cond_g <- conditional(spmod_g, newdata = newexdata, samples = 3000)
  ratio_g <- apply(cond_g, 1, sd) / preds_g$se.fit
  expect_true(all(ratio_g > 0.8 & ratio_g < 1.25))
  expect_true(all(abs(apply(cond_g, 1, mean) - preds_g$fit) < 0.1))
})

test_that("conditional.splm() no longer double-counts fixed effect uncertainty (regression test)", {
  # get_conditional_new_from_base_adjust() used to add an analytic
  # H %*% cov_betahat %*% t(H) term on top of conditional draws that were
  # already built from simulated beta draws and double-counting fixed effect
  # uncertainty. That term scales with the newdata design point x0 (via H),
  # so a covariate value far outside the observed range makes it dominate the
  # (otherwise small, since this point is spatially coincident with an
  # observed location) kriging variance and turning a subtle miscalibration
  # into an unmistakable one if the bug ever returns.
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_true(50 > max(exdata$x) + 10) # confirm x = 50 is well outside the observed range

  newdata_extreme <- newexdata[1, , drop = FALSE]
  newdata_extreme$x <- 50
  newdata_extreme$xcoord <- exdata$xcoord[1]
  newdata_extreme$ycoord <- exdata$ycoord[1]

  preds <- predict(spmod, newdata = newdata_extreme, se.fit = TRUE)
  set.seed(1)
  cond <- conditional(spmod, newdata = newdata_extreme, samples = 3000)
  ratio <- sd(cond) / preds$se.fit
  # the (fixed) implementation should land close to 1; the removed
  # double-counting term would have pushed this well above 1.15 here
  expect_true(ratio > 0.85 && ratio < 1.15)
})

test_that("conditional.spglm() draws are reproducible within the joint sampler", {
  set.seed(1)
  spmod_g <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  set.seed(1)
  cond_g <- conditional(spmod_g, newdata = newexdata[1:2, ], samples = 4)

  set.seed(1)
  expect_identical(conditional(spmod_g, newdata = newexdata[1:2, ], samples = 4), cond_g)
})

test_that("conditional.splm() draws are unchanged from a known-good seeded snapshot (regression test)", {
  set.seed(1)
  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  set.seed(1)
  cond <- conditional(spmod, newdata = newexdata[1:2, ], samples = 4)

  expect_equal(
    round(as.vector(cond), 4),
    c(1.0779, -1.2654, 1.8559, -0.858, 0.1271, -2.6709, 1.5338, -1.0532)
  )
})

test_that("Vecchia conditional simulation uses the stabilized target variance", {
  dat <- exdata[seq_len(12), , drop = FALSE]
  de <- 2
  fit <- splm(y ~ x, dat,
    xcoord = xcoord, ycoord = ycoord, ddf = "asymptotic",
    spcov_initial = spcov_initial("exponential",
      de = de, ie = 0, range = 1, known = "given"
    )
  )
  target <- newexdata[1, , drop = FALSE]
  samples <- 4
  V <- as.matrix(covmatrix(fit))
  C <- as.numeric(covmatrix(fit, target, cov_type = "pred.obs"))
  conditional_var <- de + 1e-4 * de - as.numeric(C %*% solve(V, C))

  set.seed(1)
  z <- matrix(rnorm(samples), nrow = 1)
  set.seed(1)
  observed <- get_conditional_vecchia(
    fit, target, matrix(0, nrow(dat), samples),
    local_list = list(order = 1L, method = "all", size = Inf),
    samples = samples
  )
  expect_equal(observed, sqrt(conditional_var) * z, tolerance = 1e-12)
})

test_that("conditional() local$approximation = 'vecchia' works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  cond1 <- conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", size = 10), samples = 50)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))

  # All-neighbor requests share the exact implementation and random draws.
  R <- 20
  set.seed(1)
  cond_exact <- conditional(spmod, newdata = newexdata, local = FALSE, samples = R)
  set.seed(1)
  cond_vecchia_all <- conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", method = "all"), samples = R)
  expect_identical(cond_exact, cond_vecchia_all)

  # neighbor-selection rules and distance/covariance truncation both run
  expect_vector(conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", size = 5, method = "distance"), samples = 20)[, 1])
  expect_vector(conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", size = 5, method = "covariance"), samples = 20)[, 1])

  # random effects/partition factor supported via covmatrix() reuse
  exdata_re <- exdata
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  newexdata_re <- newexdata
  newexdata_re$grp <- factor(sample(letters[1:4], NROW(newexdata_re), replace = TRUE))
  spmod_re <- splm(y ~ x, exdata_re, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~grp)
  cond_re <- conditional(spmod_re, newdata = newexdata_re, local = list(approximation = "vecchia", size = 10), samples = 30)
  expect_true(all(is.finite(cond_re)))

  # invalid local$approximation errors informatively
  expect_error(conditional(spmod, newdata = newexdata, local = list(approximation = "bogus")), "local\\$approximation must be")
})

test_that("conditional() local$approximation = 'vecchia' works for spglm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  exdata_pois <- exdata
  exdata_pois$count <- round(abs(exdata_pois$y) * 3)
  spmod <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  cond1 <- conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", size = 10), samples = 50)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))

  # All-neighbor requests share the exact implementation and random draws.
  R <- 20
  set.seed(1)
  cond_exact <- conditional(spmod, newdata = newexdata, local = FALSE, samples = R)
  set.seed(1)
  cond_vecchia_all <- conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", method = "all"), samples = R)
  expect_identical(cond_exact, cond_vecchia_all)

  # neighbor-selection rules and distance/covariance truncation both run
  expect_vector(conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", size = 5, method = "distance"), samples = 20)[, 1])
  expect_vector(conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", size = 5, method = "covariance"), samples = 20)[, 1])

  # type = "response"/"new" both work on top of the vecchia link-scale draws
  cond_response <- conditional(spmod, newdata = newexdata, type = "response", local = list(approximation = "vecchia", size = 10), samples = 20)
  expect_true(all(cond_response >= 0))
  cond_new <- conditional(spmod, newdata = newexdata, type = "new", local = list(approximation = "vecchia", size = 10), samples = 20)
  expect_true(all(cond_new == round(cond_new)))

  # random effects/partition factor supported via covmatrix() reuse
  exdata_re <- exdata_pois
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  newexdata_re <- newexdata
  newexdata_re$grp <- factor(sample(letters[1:4], NROW(newexdata_re), replace = TRUE))
  spmod_re <- spglm(count ~ x, exdata_re, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~grp)
  cond_re <- conditional(spmod_re, newdata = newexdata_re, local = list(approximation = "vecchia", size = 10), samples = 30)
  expect_true(all(is.finite(cond_re)))

  # invalid local$approximation errors informatively
  expect_error(conditional(spmod, newdata = newexdata, local = list(approximation = "bogus")), "local\\$approximation must be")

  # Only exact simulation retains the dense observed latent factorization.
  spmod_fake <- spmod
  spmod_fake$n <- 20000
  expect_no_message(
    conditional(spmod_fake, newdata = newexdata, local = list(approximation = "vecchia", size = 10), samples = 5)
  )
  expect_message(
    conditional(spmod_fake, newdata = newexdata, local = FALSE, samples = 5),
    "dense observed latent precision factorization"
  )
})

test_that("conditional() simulate_covparams = TRUE works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_false(is.null(spmod$vcov$cov)) # ddf = "satterthwaite" default for n <= 500

  cond1 <- conditional(spmod, newdata = newexdata, samples = 50, simulate_covparams = TRUE)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))

  cond_all <- conditional(spmod, newdata = newexdata, output = "all", samples = 30, simulate_covparams = TRUE)
  expect_named(cond_all, c("newdata", "beta", "object"))
  expect_equal(dim(cond_all$newdata), c(NROW(newexdata), 30))
  expect_equal(dim(cond_all$beta), c(length(coef(spmod)), 30))
  expect_equal(dim(cond_all$object), c(spmod$n, 30))

  # samples defaults to 500 (not 10,000) under simulate_covparams = TRUE
  expect_equal(ncol(conditional(spmod, newdata = newexdata, simulate_covparams = TRUE)), 1000)
  expect_equal(ncol(conditional(spmod, newdata = newexdata)), 1000)

  # propagating covariance parameter uncertainty should not shrink the
  # marginal variance of the conditional draws relative to holding covariance
  # parameters fixed
  R <- 6000
  set.seed(1)
  cond_cp <- conditional(spmod, newdata = newexdata, samples = R, simulate_covparams = TRUE)
  set.seed(1)
  cond_fixed <- conditional(spmod, newdata = newexdata, samples = R, simulate_covparams = FALSE)
  expect_true(all(apply(cond_cp, 1, sd) > apply(cond_fixed, 1, sd)))

  # random effects and anisotropy both work through the per-draw covmatrix()
  # reuse (object_b's coefficients substituted in)
  exdata_re <- exdata
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  newexdata_re <- newexdata
  newexdata_re$grp <- factor(sample(letters[1:4], NROW(newexdata_re), replace = TRUE))
  spmod_re <- splm(y ~ x, exdata_re, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~grp)
  cond_re <- conditional(spmod_re, newdata = newexdata_re, samples = 20, simulate_covparams = TRUE)
  expect_true(all(is.finite(cond_re)))

  spmod_anis <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, anisotropy = TRUE)
  cond_anis <- conditional(spmod_anis, newdata = newexdata, samples = 20, simulate_covparams = TRUE)
  expect_true(all(is.finite(cond_anis)))

  # forced back to FALSE (with a message) whenever a big-data local
  # approximation is actually active
  expect_message(
    conditional(spmod, newdata = newexdata, samples = 20, simulate_covparams = TRUE,
      local = list(approximation = "low-rank", method_base = "base", size_base = 40)
    ),
    "simulate_covparams = TRUE is not used"
  )
  # local = FALSE resolves to the exact path and leaves simulate_covparams alone
  expect_no_message(
    conditional(spmod, newdata = newexdata, samples = 20, simulate_covparams = TRUE, local = FALSE)
  )

  # a strongly-worded warning is issued above n = 500 (faking n avoids fitting
  # an actually-large model just for this check, matching the spmod_fake$n <-
  # 20000 pattern used above for vecchia's var_adj message)
  spmod_fake <- spmod
  spmod_fake$n <- 600
  expect_warning(
    conditional(spmod_fake, newdata = newexdata, samples = 2, simulate_covparams = TRUE),
    "exceedingly long"
  )
  expect_no_warning(
    conditional(spmod, newdata = newexdata, samples = 2, simulate_covparams = TRUE)
  )

  # requires object$vcov$cov (i.e. ddf = "satterthwaite") -- a clear error,
  # not a cryptic one, when that covariance matrix was never computed
  spmod_asymp <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, ddf = "asymptotic")
  expect_error(
    conditional(spmod_asymp, newdata = newexdata, samples = 2, simulate_covparams = TRUE),
    "requires object\\$vcov\\$cov"
  )

  # simulate_covparams must be a single logical
  expect_error(conditional(spmod, newdata = newexdata, simulate_covparams = "yes"), "simulate_covparams must be")
})

test_that("conditional() output = 'cov'/'spcov'/'randcov' works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  spcov_names <- names(coef(spmod, type = "spcov"))

  cov_out <- conditional(spmod, newdata = newexdata, samples = 30, simulate_covparams = TRUE, output = "cov")
  expect_equal(dim(cov_out), c(length(spcov_names), 30))
  expect_equal(rownames(cov_out), spcov_names)
  expect_true(all(is.finite(cov_out)))

  # no random effects in this model -- "spcov" and "cov" coincide, "randcov" is NULL
  spcov_out <- conditional(spmod, newdata = newexdata, samples = 30, simulate_covparams = TRUE, output = "spcov")
  expect_equal(dim(spcov_out), c(length(spcov_names), 30))
  randcov_out <- conditional(spmod, newdata = newexdata, samples = 30, simulate_covparams = TRUE, output = "randcov")
  expect_null(randcov_out)

  # combining with the pre-existing outputs still works, in any combination
  combo <- conditional(spmod, newdata = newexdata,
    samples = 30, simulate_covparams = TRUE, output = c("newdata", "beta", "cov", "spcov")
  )
  expect_named(combo, c("newdata", "beta", "cov", "spcov"))
  expect_equal(dim(combo$newdata), c(NROW(newexdata), 30))
  expect_equal(dim(combo$beta), c(length(coef(spmod)), 30))
  expect_equal(dim(combo$cov), c(length(spcov_names), 30))
  expect_equal(dim(combo$spcov), c(length(spcov_names), 30))

  # requesting these without simulate_covparams = TRUE errors clearly
  expect_error(
    conditional(spmod, newdata = newexdata, samples = 5, output = "cov"),
    "only include \"cov\", \"spcov\", or \"randcov\" when simulate_covparams = TRUE"
  )

  # random effects: "randcov" is populated, and "cov" stacks spcov then randcov rows
  exdata_re <- exdata
  exdata_re$grp <- factor(sample(letters[1:4], NROW(exdata_re), replace = TRUE))
  newexdata_re <- newexdata
  newexdata_re$grp <- factor(sample(letters[1:4], NROW(newexdata_re), replace = TRUE))
  spmod_re <- splm(y ~ x, exdata_re, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~grp)
  randcov_names <- names(coef(spmod_re, type = "randcov"))

  randcov_out2 <- conditional(spmod_re, newdata = newexdata_re, samples = 30, simulate_covparams = TRUE, output = "randcov")
  expect_equal(dim(randcov_out2), c(length(randcov_names), 30))
  expect_equal(rownames(randcov_out2), randcov_names)

  cov_out2 <- conditional(spmod_re, newdata = newexdata_re, samples = 30, simulate_covparams = TRUE, output = "cov")
  expect_equal(rownames(cov_out2), c(names(coef(spmod_re, type = "spcov")), randcov_names))
})

test_that("simulate_theta_draw()'s clamp-to-boundary fallback is always valid", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  spcov_params_full <- coef(spmod, type = "spcov")
  spcov_type <- class(spcov_params_full)
  spcov_names_free <- names(spmod$is_known$spcov)[!spmod$is_known$spcov]
  theta_hat_free <- spcov_params_full[spcov_names_free]

  # a huge covariance matrix forces essentially every reject-and-redraw
  # attempt to fail, exercising the clamp-to-boundary fallback
  huge_vcov <- diag(1e6, length(theta_hat_free))
  dimnames(huge_vcov) <- list(spcov_names_free, spcov_names_free)
  huge_lowchol <- t(chol(huge_vcov))

  set.seed(1)
  draws <- lapply(1:30, function(i) {
    simulate_theta_draw(
      theta_hat_free, huge_lowchol, spcov_type, spcov_names_free,
      character(0), spcov_params_full, NULL,
      max_attempts = 5
    )
  })
  de_vals <- vapply(draws, function(x) x$spcov_params_b[["de"]], numeric(1))
  ie_vals <- vapply(draws, function(x) x$spcov_params_b[["ie"]], numeric(1))
  range_vals <- vapply(draws, function(x) x$spcov_params_b[["range"]], numeric(1))

  expect_true(all(is.finite(c(de_vals, ie_vals, range_vals))))
  expect_true(all(de_vals >= 0) && all(ie_vals >= 0) && all(range_vals >= 0))
  # a clamped range of exactly 0 must also force de to 0 (avoids
  # division-by-zero in exp(-d/range)-style formulas downstream)
  expect_true(all(de_vals[range_vals == 0] == 0))
})



# Additional checks
test_that("conditional works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  cond1 <- conditional(spmod, newdata = newexdata, samples = 50)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))

  cond_all <- conditional(spmod, newdata = newexdata, output = "all", samples = 50)
  expect_type(cond_all, "list")
  expect_named(cond_all, c("newdata", "beta", "object"))
  expect_equal(dim(cond_all$newdata), c(NROW(newexdata), 50))
  expect_equal(dim(cond_all$beta), c(length(coef(spmod)), 50))
  expect_equal(dim(cond_all$object), c(spmod$n, 50))

  # falls back to object$newdata (the missing rows) when newdata is omitted
  exdata_miss <- exdata
  exdata_miss$y[1:5] <- NA
  spmod_miss <- splm(y ~ x, exdata_miss, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  cond_miss <- conditional(spmod_miss, samples = 50)
  expect_equal(nrow(cond_miss), 5)
})

test_that("conditional works for spglm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  exdata_pois <- exdata
  exdata_pois$count <- round(abs(exdata_pois$y) * 3)
  spmod <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  cond_link <- conditional(spmod, newdata = newexdata, samples = 50)
  expect_equal(dim(cond_link), c(NROW(newexdata), 50))

  cond_response <- conditional(spmod, newdata = newexdata, type = "response", samples = 50)
  expect_true(all(cond_response >= 0))

  cond_new <- conditional(spmod, newdata = newexdata, type = "new", samples = 50)
  expect_true(all(cond_new >= 0))
  expect_equal(cond_new, round(cond_new)) # poisson draws are counts

  for (samples in c(1, 5)) {
    for (type in c("link", "response", "new")) {
      draws <- conditional(spmod, newdata = newexdata[1, , drop = FALSE],
        type = type, samples = samples)
      expect_equal(dim(draws), c(1, samples))
      expect_true(all(is.finite(draws)))
      if (type == "new") {
        expect_true(all(draws >= 0 & draws == floor(draws)))
      }
    }
  }
})

test_that("conditional() local$approximation = 'vecchia' works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  cond1 <- conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", size = 10), samples = 50)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))
})

test_that("conditional() local$approximation = 'vecchia' works for spglm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  exdata_pois <- exdata
  exdata_pois$count <- round(abs(exdata_pois$y) * 3)
  spmod <- spglm(count ~ x, exdata_pois, family = "poisson", spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)

  cond1 <- conditional(spmod, newdata = newexdata, local = list(approximation = "vecchia", size = 10), samples = 50)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))
})

test_that("conditional() simulate_covparams = TRUE works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_false(is.null(spmod$vcov$cov)) # ddf = "satterthwaite" default for n <= 500

  cond1 <- conditional(spmod, newdata = newexdata, samples = 50, simulate_covparams = TRUE)
  expect_true(is.matrix(cond1))
  expect_equal(dim(cond1), c(NROW(newexdata), 50))
  expect_true(all(is.finite(cond1)))
})

test_that("conditional() output = 'cov'/'spcov'/'randcov' works for splm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod <- splm(y ~ x, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  spcov_names <- names(coef(spmod, type = "spcov"))

  cov_out <- conditional(spmod, newdata = newexdata, samples = 30, simulate_covparams = TRUE, output = "cov")
  expect_equal(dim(cov_out), c(length(spcov_names), 30))
  expect_equal(rownames(cov_out), spcov_names)
  expect_true(all(is.finite(cov_out)))
})


# Additional areal checks

test_that("areal simulation uses the fitted missing rows and output conventions", {
  for (family in list(NULL, "poisson")) {
    fit <- areal_conditional_fixture(family)
    set.seed(1)
    draws <- conditional(fit, output = "all", samples = 8)
    set.seed(1)
    expect_identical(conditional(fit, newdata = fit$newdata, output = "all", samples = 8), draws)
    expect_named(draws, c("newdata", "beta", "object"))
    expect_equal(dim(draws$newdata), c(3, 8))
    expect_equal(dim(draws$beta), c(2, 8))
    expect_identical(rownames(draws$newdata), as.character(fit$missing_index))
    expect_identical(rownames(draws$object), as.character(fit$observed_index))
    expect_identical(rownames(draws$beta), names(coef(fit)))
    observed <- if (is.null(family)) model.response(model.frame(fit)) else fitted(fit, type = "link")
    expect_equal(unname(draws$object), matrix(rep(observed, 8), ncol = 8))
    expect_identical(conditional(fit, output = "object", samples = 8), draws$object)
    for (output in list("newdata", "beta", c("object", "beta"))) {
      set.seed(1)
      expected <- if (length(output) == 1L) draws[[output]] else draws[output]
      expect_identical(conditional(fit, output = output, samples = 8), expected)
    }
    expect_error(conditional(fit, newdata = fit$newdata[3:1, ]), "newdata cannot")
    expect_error(conditional(fit, newdata = NULL), "newdata cannot")
    for (local in list(FALSE, TRUE, list(approximation = "vecchia"))) {
      expect_error(conditional(fit, local = local), "only exact")
    }
    expect_error(conditional(fit, simulate_covparams = TRUE), "not supported")
    expect_error(conditional(fit, output = "cov"), "output must")
    expect_error(conditional(fit, output = character()), "output must")
    for (samples in list(0, -1, 1.5, NA, Inf, c(1, 2))) {
      expect_error(conditional(fit, samples = samples), "positive integer")
    }
    complete <- areal_conditional_fixture(family, missing = integer())
    expect_error(conditional(complete), "No missing data")
    expect_error(conditional(complete, newdata = NULL), "No missing data")
  }
})

test_that("areal simulation preserves single-row and single-draw dimensions", {
  for (family in list(NULL, "poisson")) {
    fit <- areal_conditional_fixture(family, spcov_type = "sar", missing = 7,
      formula = y ~ 1 + offset(off))
    for (samples in c(1, 5)) {
      draws <- conditional(fit, samples = samples, output = "all")
      expect_equal(dim(draws$newdata), c(1, samples))
      expect_equal(dim(draws$beta), c(1, samples))
      if (!is.null(family)) {
        for (type in c("response", "new")) {
          expect_equal(dim(conditional(fit, samples = samples, type = type)), c(1, samples))
        }
      }
    }
  }
})

test_that("areal joint preparation agrees with an explicitly supplied covariance", {
  for (spcov_type in c("car", "sar")) {
    fit <- areal_conditional_fixture("poisson", spcov_type)
    expect_equal(get_conditional_glm_joint(fit),
      get_conditional_glm_joint(fit, cov_lowchol = t(chol(covmatrix(fit)))))
  }
})

test_that("areal GLM transformations and binomial sizes match link draws", {
  fit <- areal_conditional_fixture("poisson")
  set.seed(1)
  link <- conditional(fit, samples = 20)
  set.seed(1)
  expect_equal(conditional(fit, samples = 20, type = "response"), exp(link))
  set.seed(1)
  expect_equal(conditional(fit, fit$newdata, "newdata", "response", 20), exp(link))
  counts <- conditional(fit, samples = 20, type = "new")
  expect_true(all(counts >= 0 & counts == floor(counts)))

  fit <- areal_conditional_fixture("binomial")
  expect_identical(unname(fit$missing_index), c(2L, 7L, 15L))
  size <- c(1, 5, 9)
  set.seed(1)
  link <- conditional(fit, samples = 20)
  set.seed(1)
  expect_equal(conditional(fit, samples = 20, type = "response", newdata_size = size),
    sweep(plogis(link), 1, size, "*"))
  counts <- conditional(fit, samples = 20, type = "new", newdata_size = size)
  expect_true(all(counts >= 0 & counts == floor(counts) & counts <= size))
  set.seed(1)
  scalar <- conditional(fit, samples = 3, type = "response", newdata_size = 5)
  set.seed(1)
  expect_equal(conditional(fit, samples = 3, type = "response", newdata_size = rep(5, 3)), scalar)
  for (size in list(c(1, 2), -1, 1.5, NA, Inf)) {
    expect_error(conditional(fit, newdata_size = size), "newdata_size must")
  }
})

test_that("areal simulation reuses preparation for all sample columns", {
  prepare <- get_conditional_areal_cov
  counts <- 0L
  local_mocked_bindings(get_conditional_areal_cov = function(object) {
    counts <<- counts + 1L
    prepare(object)
  })
  for (family in list(NULL, "poisson")) {
    fit <- areal_conditional_fixture(family)
    for (samples in c(1, 100)) {
      counts <- 0L
      conditional(fit, samples = samples)
      expect_identical(counts, 1L)
    }
    counts <- 0L
    conditional(fit, samples = 100, output = "object")
    expect_identical(counts, 0L)
  }
})


# Exact GLM preparation and factor reuse

test_that("exact GLM preparation uses corrected coefficients and two dense factors", {
  fixture <- joint_fixture()
  fit <- fixture$fit
  sizes <- integer()
  local_mocked_bindings(chol = function(x, ...) {
    sizes <<- c(sizes, NROW(x))
    base::chol(x, ...)
  })
  joint <- get_conditional_glm_joint(fit)
  expect_identical(sizes, c(18L, 18L, 2L))
  expect_equal(tcrossprod(joint$cov_betahat_lowchol), vcov(fit))
  sizes <- integer()
  for (samples in c(1, 100)) draw_conditional_glm_joint(joint, samples)
  expect_length(sizes, 0L)
  fit$local_index <- rep(1:3, 6)
  expect_identical(get_conditional_glm_joint(fit), joint)
  sizes <- integer()
  expect_equal(get_conditional_glm_joint(fit, joint$cov_lowchol), joint)
  expect_identical(sizes, c(18L, 2L))
})

test_that("complete exact GLM calls factor observed covariance only once", {
  point <- joint_fixture()
  areal <- areal_conditional_fixture("poisson")
  original_draw <- draw_conditional_glm_joint
  sizes <- integer()
  local_mocked_bindings(chol = function(x, ...) {
    sizes <<- c(sizes, NROW(x))
    base::chol(x, ...)
  }, draw_conditional_glm_joint = function(joint, samples, residual = FALSE) {
    if (!residual) stop("exact callers need residual draws directly")
    original_draw(joint, samples, residual = residual)
  })
  for (samples in c(1, 100)) {
    sizes <- integer()
    conditional(point$fit, point$newdata, local = FALSE, samples = samples)
    expect_identical(sizes, c(18L, 18L, 2L, 4L))
    sizes <- integer()
    conditional(areal, samples = samples)
    # CAR covariance construction also factors the full 18-site graph precision.
    expect_identical(sizes, c(18L, 15L, 3L, 15L, 2L))
  }
})

test_that("direct residual draws retain the latent and coefficient joint", {
  fixture <- joint_fixture()
  joint <- get_conditional_glm_joint(fixture$fit)
  for (samples in c(1, 7)) {
    set.seed(1)
    latent <- draw_conditional_glm_joint(joint, samples)
    set.seed(1)
    residual <- draw_conditional_glm_joint(joint, samples, residual = TRUE)
    expect_identical(residual$beta, latent$beta)
    expect_equal(unname(residual$residual),
      unname(latent$w - model.matrix(fixture$fit) %*% latent$beta), tolerance = 1e-12)
  }
})


# Low-rank GLM preparation

test_that("low-rank GLM simulation avoids full observed preparation and fitting blocks", {
  fixture <- joint_fixture()
  fit <- fixture$fit
  local <- list(approximation = "low-rank", size_base = 9, reorder_base = "none", size_new = 2,
    reorder_new = "none", chunk_size = 3)
  local_mocked_bindings(get_conditional_glm_joint = function(...) stop("dense preparation"))
  covariance <- covmatrix.spglm
  dimensions <- integer()
  local_mocked_bindings(covmatrix.spglm = function(object, ...) {
    dimensions <<- c(dimensions, NROW(object$obdata))
    covariance(object, ...)
  })
  set.seed(1)
  draws <- conditional(fit, fixture$newdata, samples = 7, local = local, output = "all")
  expect_lte(max(dimensions), 9)
  expect_equal(dim(draws$newdata), c(4, 7))
  expect_equal(dim(draws$beta), c(2, 7))
  expect_equal(draws$object[, 1], unname(fitted(fit, type = "link")))
  fit$local_index <- rep(1:3, 6)
  set.seed(1)
  expect_identical(conditional(fit, fixture$newdata, samples = 7, local = local, output = "all"), draws)
})

test_that("low-rank GLM factors are reused across draw counts and blocks", {
  fixture <- joint_fixture()
  prepare_base <- get_conditional_glm_base
  prepare_block <- get_conditional_glm_block
  counts <- c(base = 0L, block = 0L)
  local_mocked_bindings(
    get_conditional_glm_base = function(...) {
      counts[["base"]] <<- counts[["base"]] + 1L
      prepare_base(...)
    },
    get_conditional_glm_block = function(...) {
      counts[["block"]] <<- counts[["block"]] + 1L
      prepare_block(...)
    })
  for (samples in c(1, 100)) {
    for (chunk_size in c(1, 1000)) {
      counts[] <- 0L
      conditional(fixture$fit, fixture$newdata, samples = samples,
        local = list(approximation = "low-rank", size_base = 9, reorder_base = "none", size_new = 2,
          reorder_new = "none", chunk_size = chunk_size))
      expect_identical(counts, c(base = 1L, block = 2L))
    }
  }
})

test_that("rank-deficient bases work and invalid chunk sizes error", {
  fixture <- joint_fixture()
  local <- list(approximation = "low-rank", size_base = 1, reorder_base = "none")
  expect_true(all(is.finite(conditional(fixture$fit, fixture$newdata, samples = 2, local = local))))
  local$size_base <- 9
  for (size in list(0, NA, Inf, 1.5, c(1, 2))) {
    local$chunk_size <- size
    expect_error(conditional(fixture$fit, fixture$newdata, local = local), "chunk_size must")
  }
})

test_that("full-base joint distribution matches the global exact approximation", {
  fixture <- joint_fixture()
  fit <- fixture$fit
  base <- get_conditional_glm_base(fit, seq_len(fit$n))
  joint <- get_conditional_glm_joint(fit)
  H <- base$cov_lowchol_base %*% chol2inv(t(base$cov_lowchol_mH)) %*% t(base$cov_lowchol_base)
  M <- base$cov_lowchol_base %*% base$wts_latent
  H_exact <- chol2inv(joint$cond_prec_upchol)
  M_exact <- joint$wts_latent
  expect_equal(unname(H), unname(H_exact), tolerance = 1e-7)
  expect_equal(unname(M), unname(M_exact), tolerance = 1e-7)
  target <- lowrank_joint_reference(fit, fixture$newdata, seq_len(fit$n))
  exact <- joint_reference(fit, fixture$newdata)
  expect_equal(target$mean, exact$mean, tolerance = 1e-7)
  expect_equal(unname(target$covariance), unname(exact$covariance), tolerance = 1e-7)
})

test_that("low-rank batching supports kmeans blocks and near-singular covariance", {
  for (near in c(FALSE, TRUE)) {
    fixture <- joint_fixture(near = near)
    draws <- conditional(fixture$fit, fixture$newdata, samples = 5,
      local = list(approximation = "low-rank", size_base = 9, reorder_base = "none", size_new = 2,
        kmeans_new = TRUE, chunk_size = 2))
    expect_equal(dim(draws), c(4, 5))
    expect_true(all(is.finite(draws)))
  }
})

test_that("low-rank retains coefficient uncertainty when the base misses a factor level", {
  fixture <- joint_fixture()
  fit <- spglm(y ~ group + offset(off), fixture$fit$obdata, family = "poisson",
    xcoord = "cx", ycoord = "cy", contrasts = list(group = "contr.sum"),
    spcov_initial = spcov_initial("exponential", de = 0.4, ie = 0.2, range = 0.3, known = "given"))
  base <- get_conditional_glm_base(fit, 1L)
  expect_equal(unname(tcrossprod(base$cov_betahat_lowchol)), unname(vcov(fit)))
  target <- lowrank_joint_reference(fit, fixture$newdata, 1L)
  expect_equal(unname(target$covariance[1:3, 1:3]), unname(vcov(fit)))
  expect_true(all(is.finite(conditional(fit, fixture$newdata, samples = 2,
    local = list(approximation = "low-rank", size_base = 1, reorder_base = "none")))))
})

test_that("low-rank releases prediction factors before preparing later blocks", {
  fixture <- joint_fixture()
  prepare <- get_conditional_glm_block
  events <- character()
  local_mocked_bindings(get_conditional_glm_block = function(...) {
    events <<- c(events, "prepare")
    prepare(...)
  }, rnorm = function(n, ...) {
    events <<- c(events, "draw")
    stats::rnorm(n, ...)
  })
  conditional(fixture$fit, fixture$newdata, samples = 5,
    local = list(approximation = "low-rank", size_base = 9, reorder_base = "none", size_new = 1,
      reorder_new = "none", chunk_size = 2))
  preparation <- which(events == "prepare")
  expect_length(preparation, 4)
  expect_equal(diff(preparation), rep(4L, 3))
})

test_that("resolved conditional settings determine exact dispatch", {
  fixture <- joint_fixture()
  for (setting in list(NULL, FALSE)) {
    expect_true(get_local_list_conditional(setting, fixture$fit, fixture$newdata)$exact)
  }
  expect_false(get_local_list_conditional(list(approximation = "low-rank", method_base = "all", method_new = "all"),
    fixture$fit, fixture$newdata)$exact)
  fixture$fit$n <- 5001
  local_mocked_bindings(get_local_list_conditional_lowrank = function(local, ...) local)
  expect_message(local <- get_local_list_conditional(NULL, fixture$fit, fixture$newdata), "local = TRUE")
  expect_false(local$exact)
})

test_that("parallel low-rank preparation bounds the number of blocks", {
  fixture <- joint_fixture()
  single <- get_local_list_conditional(list(approximation = "low-rank", method_new = "all", parallel = TRUE, ncores = 2),
    fixture$fit, fixture$newdata)
  expect_equal(single$ncores, 1)
  payloads <- integer()
  local_mocked_bindings(.package = "parallel", makeCluster = function(...) NULL,
    stopCluster = function(...) NULL,
    parLapply = function(cl, X, fun, ...) {
      payloads <<- c(payloads, length(X))
      lapply(X, fun, ...)
    })
  conditional(fixture$fit, fixture$newdata, samples = 3,
    local = list(approximation = "low-rank", size_base = 9, reorder_base = "none", size_new = 1,
      reorder_new = "none", parallel = TRUE, ncores = 2, chunk_size = 1))
  expect_identical(payloads, c(2L, 2L))
})


# Revised observed-latent preparation is tested in test-extras-conditional-vecchia.R.

# Areal joint moments

areal_conditional_reference <- function(fit) {
  X <- model.matrix(fit)
  mf <- model.frame(delete.response(terms(fit)), fit$newdata, xlev = fit$xlevels)
  Xnew <- model.matrix(delete.response(terms(fit)), mf, contrasts.arg = fit$contrasts)
  Xnew <- Xnew[, colnames(X), drop = FALSE]
  offset <- model.offset(mf)
  P <- solve(covmatrix(fit))
  C <- covmatrix(fit, fit$newdata)
  A <- C %*% P
  R <- covmatrix(fit, fit$newdata, cov_type = "pred.pred") - A %*% t(C)
  J <- Xnew - A %*% X
  if (inherits(fit, "spgautor")) {
    eta <- fitted(fit, type = "link")
    w <- eta - model.offset(model.frame(fit))
    curvature <- switch(fit$family,
      poisson = exp(eta),
      binomial = rowSums(model.response(model.frame(fit))) * plogis(eta) * (1 - plogis(eta)))
    H <- solve(P + diag(curvature))
    M <- H %*% P %*% X
    F <- J + A %*% M
    beta_cov <- vcov(fit)
    pred_cov <- R + A %*% H %*% t(A) + F %*% beta_cov %*% t(F)
    cross <- beta_cov %*% t(F)
  } else {
    w <- model.response(model.frame(fit)) - model.offset(model.frame(fit))
    beta_cov <- vcov(fit)
    pred_cov <- R + J %*% beta_cov %*% t(J)
    cross <- beta_cov %*% t(J)
  }
  mu <- Xnew %*% coef(fit) + A %*% (w - X %*% coef(fit)) + offset
  list(mean = c(coef(fit), mu),
    covariance = rbind(cbind(beta_cov, cross), cbind(t(cross), pred_cov)))
}

check_areal_conditional_moments <- function(fit) {
  target <- areal_conditional_reference(fit)
  p <- length(coef(fit))
  prediction <- predict(fit, se.fit = TRUE)
  expect_equal(as.numeric(prediction$fit), unname(target$mean[-seq_len(p)]), tolerance = 1e-8)
  # Stored vcov uses the final fitting Hessian, which can precede the last
  # latent update; the simulation evaluates conditional curvature at fitted w.
  expect_equal(as.numeric(prediction$se.fit)^2,
    unname(diag(target$covariance)[-seq_len(p)]), tolerance = 1e-6)
  expect_equal(unname(vcov(fit)), unname(target$covariance[seq_len(p), seq_len(p), drop = FALSE]), tolerance = 1e-7)
  set.seed(1)
  draws <- conditional(fit, output = c("beta", "newdata"), samples = 4000)
  draws <- rbind(draws$beta, draws$newdata)
  S <- target$covariance
  mcse <- sqrt((S^2 + outer(diag(S), diag(S))) / (NCOL(draws) - 1))
  expect_lt(max(abs(cov(t(draws)) - S) / mcse), 6)
  expect_lt(max(abs(rowMeans(draws) - target$mean) / sqrt(diag(S) / NCOL(draws))), 6)
}

test_that("CAR and SAR conditional moments agree with analytical prediction", {
  for (family in list(NULL, "poisson")) {
    for (spcov_type in c("car", "sar")) {
      for (row_st in c(FALSE, TRUE)) {
        fit <- areal_conditional_fixture(family, spcov_type, row_st)
        check_areal_conditional_moments(fit)
      }
    }
  }
})

test_that("areal joint simulation retains random slopes, partitions and contrasts", {
  for (family in list(NULL, "poisson")) {
    for (partition in c(FALSE, TRUE)) {
      fit <- areal_conditional_fixture(family, random = TRUE, partition = partition,
        formula = y ~ x + group + offset(off), polygons = TRUE)
      check_areal_conditional_moments(fit)
      expect_identical(rownames(conditional(fit, output = "beta", samples = 2)), names(coef(fit)))
    }
  }
})

test_that("binomial areal joint moments account for trial sizes greater than one", {
  for (spcov_type in c("car", "sar")) {
    fit <- areal_conditional_fixture("binomial", spcov_type)
    expect_true(all(rowSums(model.response(model.frame(fit))) > 1))
    check_areal_conditional_moments(fit)
  }
})

test_that("areal simulation supports dispersion GLM families", {
  for (family in c("nbinomial", "Gamma", "inverse.gaussian", "beta")) {
    fit <- areal_conditional_fixture(family)
    set.seed(1)
    link <- conditional(fit, samples = 20)
    set.seed(1)
    response <- conditional(fit, type = "response", samples = 20)
    expect_equal(response, if (family == "beta") plogis(link) else exp(link))
    new <- conditional(fit, type = "new", samples = 20)
    expect_equal(dim(new), c(3, 20))
    expect_true(all(is.finite(new) & new >= 0))
    if (family == "beta") expect_true(all(new < 1))
    if (family == "nbinomial") expect_equal(new, floor(new))
  }
})


# Exact and local GLM joint moments

test_that("exact joint factors and full sampled moments match independent matrices", {
  fixture <- joint_fixture()
  fit <- fixture$fit
  target <- joint_reference(fit, fixture$newdata)
  joint <- get_conditional_glm_joint(fit)
  expect_equal(unname(joint$wts_latent), unname(target$M), tolerance = 1e-10)
  expect_equal(chol2inv(joint$cond_prec_upchol), unname(target$H), tolerance = 1e-10)
  expect_equal(unname(tcrossprod(joint$cov_betahat_lowchol)), unname(target$C), tolerance = 1e-10)
  expect_equal(unname(vcov(fit)), unname(target$covariance[1:2, 1:2]), tolerance = 1e-4)
  prediction <- predict(fit, fixture$newdata, se.fit = TRUE, type = "link", local = FALSE)
  expect_equal(as.numeric(prediction$fit), unname(target$mean[-(1:2)]), tolerance = 1e-8)
  expect_equal(as.numeric(prediction$se.fit)^2, unname(diag(target$covariance)[-(1:2)]), tolerance = 1e-8)
  excess <- target$F %*% target$C %*% t(target$F)
  expect_gt(max(abs(excess)), 0.001)
  set.seed(1)
  draws <- conditional(fit, fixture$newdata, output = "all", samples = 30000, local = FALSE)
  expect_joint_moments(rbind(draws$beta, draws$newdata), target)
  expect_equal(draws$object[, 1], unname(fitted(fit, type = "link")))
  expect_identical(rownames(draws$beta), names(coef(fit)))
})

test_that("local fitting and simulation preserve the constructed joint moments", {
  simulation <- list(FALSE,
    list(approximation = "low-rank", method_base = "all", method_new = "all"),
    list(approximation = "low-rank", method_base = "base", size_base = 9, reorder_base = "none", method_new = "base", size_new = 2, reorder_new = "none"),
    list(approximation = "vecchia", method = "all", ordering = "none"),
    list(approximation = "vecchia", method = "distance", size = 5, ordering = "none"))
  for (adjustment in c("exact", "none", "theoretical", "pooled", "empirical")) {
    local_fit <- if (adjustment == "exact") FALSE else list(index = rep(c(3, 1, 2, 1, 2, 3), 3), var_adjust = adjustment)
    fixture <- joint_fixture(local_fit)
    fit <- fixture$fit
    full <- joint_reference(fit, fixture$newdata)
    for (reference in list(
        lowrank_joint_reference(fit, fixture$newdata, seq_len(fit$n)),
        vecchia_joint_reference(fit, fixture$newdata))) {
      expect_equal(reference$mean, full$mean, tolerance = 1e-9)
      expect_equal(unname(reference$covariance), unname(full$covariance), tolerance = 1e-9)
    }
    for (i in seq_along(simulation)) {
      target <- joint_reference(fit, fixture$newdata,
        base = if (i == 3) 1:9 else 1:18,
        blocks = if (i == 3) list(1:2, 3:4) else list(1:4),
        neighbors = if (i == 5) 5 else NULL, base_latent = i %in% c(2, 3))
      if (i %in% c(4, 5)) {
        target <- vecchia_joint_reference(fit, fixture$newdata, size = if (i == 4) Inf else 5)
      }
      joint <- get_conditional_glm_joint(fit)
      exact_target <- joint_reference(fit, fixture$newdata)
      expect_equal(unname(joint$wts_latent), unname(exact_target$M), tolerance = 1e-9)
      expect_equal(chol2inv(joint$cond_prec_upchol), unname(exact_target$H), tolerance = 1e-9)
      expect_equal(unname(tcrossprod(joint$cov_betahat_lowchol)),
        unname(exact_target$C), tolerance = 1e-9)
      expect_equal(unname(vcov(fit)), unname(exact_target$covariance[1:2, 1:2]), tolerance = 1e-4)
      set.seed(1)
      draws <- conditional(fit, fixture$newdata, output = c("beta", "newdata"), samples = 30000, local = simulation[[i]])
      expect_joint_moments(rbind(draws$beta, draws$newdata), target)
    }
  }
})

test_that("exact observed latent and coefficient draws have the specified joint", {
  fixture <- joint_fixture(local = list(index = rep(1:3, 6), var_adjust = "theoretical"))
  fit <- fixture$fit
  reference <- joint_reference(fit, fixture$newdata)
  C <- reference$C
  M <- reference$M
  H <- reference$H
  target <- list(mean = c(coef(fit), fitted(fit, type = "link") - fit$obdata$off),
    covariance = rbind(cbind(C, C %*% t(M)),
      cbind(M %*% C, H + M %*% C %*% t(M))))
  set.seed(1)
  draws <- draw_conditional_glm_joint(get_conditional_glm_joint(fit), 4000)
  expect_joint_moments(rbind(draws$beta, draws$w), target)
})

test_that("coupled draws preserve family, offset, dimension", {
  for (family in c("poisson", "nbinomial", "binomial", "Gamma", "inverse.gaussian", "beta")) {
    fixture <- joint_fixture(family = family)
    for (local in list(FALSE, list(approximation = "vecchia", size = 5, ordering = "none"),
        list(approximation = "low-rank", size_base = 9, reorder_base = "none", chunk_size = 2))) {
      fit <- fixture$fit
      new <- fixture$newdata[1, , drop = FALSE]
      set.seed(1)
      link <- conditional(fit, new, samples = 1, local = local, output = "all")
      set.seed(1)
      response <- conditional(fit, new, samples = 1, local = local, type = "response", newdata_size = 5)
      expect_equal(dim(link$newdata), c(1L, 1L))
      expect_equal(dim(response), c(1L, 1L))
      expect_equal(as.numeric(response), as.numeric(invlink(link$newdata, family, size = 5)))
      expect_equal(link$object[, 1], unname(fitted(fit, type = "link")))
      set.seed(1)
      expect_identical(conditional(fit, new, samples = 1, local = local, output = "all"), link)
      set.seed(1)
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
      list(approximation = "low-rank", method_base = "base", size_base = 9, reorder_base = "none",
        method_new = "base", size_new = 2, reorder_new = "none", parallel = TRUE, ncores = 2),
      list(approximation = "vecchia", method = "distance", size = 5, ordering = "none"))) {
    lowrank <- is.list(simulation) && isTRUE(simulation$parallel)
    vecchia <- is.list(simulation) && identical(simulation$approximation, "vecchia")
    target <- joint_reference(fit, fixture$newdata, base = if (lowrank) 1:9 else 1:18,
      blocks = if (lowrank) list(1:2, 3:4) else list(1:4), neighbors = if (vecchia) 5 else NULL,
      base_latent = lowrank)
    if (vecchia) target <- vecchia_joint_reference(fit, fixture$newdata, size = 5)
    set.seed(1)
    draws <- conditional(fit, fixture$newdata, output = "all", samples = 30000, local = simulation)
    expect_joint_moments(rbind(draws$beta, draws$newdata), target)
    expect_identical(fit, original)
  }
})

test_that("coordinate forms, missing responses and intercept-only fits remain supported", {
  for (coordinates in c("data", "sf", "one")) {
    fixture <- joint_fixture(coordinates = coordinates, intercept = TRUE, missing = TRUE)
    fit <- fixture$fit
    for (local in list(FALSE, list(approximation = "vecchia", method = "all", ordering = "none"),
        list(approximation = "low-rank", size_base = 9, reorder_base = "none", chunk_size = 2))) {
      set.seed(1)
      explicit <- conditional(fit, fit$newdata, output = "all", samples = 1, local = local)
      set.seed(1)
      implicit <- conditional(fit, output = "all", samples = 1, local = local)
      expect_identical(explicit, implicit)
      expect_equal(dim(explicit$beta), c(1L, 1L))
      expect_equal(dim(explicit$newdata), c(2L, 1L))
      expect_equal(explicit$object[, 1], unname(fitted(fit, type = "link")))
      expect_true(all(is.finite(explicit$newdata)))
    }
  }
})

test_that("near-singular spatial covariance and invalid latent precision are properly handled", {
  fixture <- joint_fixture(near = TRUE)
  expect_true(all(is.finite(conditional(fixture$fit, fixture$newdata, samples = 20))))
  fixture <- joint_fixture()
  testthat::local_mocked_bindings(get_D = function(...) diag(1e6, 18))
  expect_error(get_conditional_glm_joint(fixture$fit), "latent precision is not positive definite")
})

test_that("factor contrasts and parallel draws retain their output structure", {
  fixture <- joint_fixture()
  data <- fixture$fit$obdata
  data$x <- seq(-1, 1, length.out = NROW(data))
  fit <- spglm(y ~ group + x + offset(off), family = "poisson", data = data,
    xcoord = "cx", ycoord = "cy", contrasts = list(group = "contr.sum"),
    spcov_initial = spcov_initial("exponential", de = 0.4, ie = 0.2, range = 0.3, known = "given"))
  newdata <- fixture$newdata[c(3, 2, 1, 4), ]
  set.seed(1)
  draws <- conditional(fit, newdata, output = "all", samples = 30000, local = FALSE)
  target <- joint_reference(fit, newdata)
  expect_joint_moments(rbind(draws$beta, draws$newdata), target)
  expect_identical(rownames(draws$beta), names(coef(fit)))
  local <- list(approximation = "low-rank", method_base = "base", size_base = 9, reorder_base = "random",
    method_new = "base", size_new = 2, reorder_new = "random", kmeans_new = FALSE,
    parallel = TRUE, ncores = 2)
  set.seed(1)
  parallel_draws <- conditional(fit, newdata, output = "all", samples = 10, local = local)
  expect_equal(dim(parallel_draws$newdata), c(NROW(newdata), 10))
  expect_true(all(is.finite(parallel_draws$newdata)))
  expect_identical(rownames(parallel_draws$beta), names(coef(fit)))
})


# Local GLM joint moments

test_that("Vecchia joint moments retain previous latent and coefficient dependence", {
  for (family in c("poisson", "binomial")) {
    fixture <- joint_fixture(family = family)
    fixture$newdata$cx <- 10 + seq_len(4) / 100
    fixture$newdata$cy <- 0
    for (size in c(1, 5)) {
      for (chunk in c(100, 1000)) {
        target <- vecchia_joint_reference(fixture$fit, fixture$newdata, size)
        set.seed(1)
        draws <- conditional(fixture$fit, fixture$newdata, samples = 4000,
          output = c("beta", "newdata"), local = list(approximation = "vecchia",
            method = "distance", size = size, ordering = "none", chunk_size = chunk))
        expect_joint_moments(rbind(draws$beta, draws$newdata), target)
      }
    }
  }
})

test_that("Vecchia covariance selection respects random slopes and partitions", {
  fixture <- joint_fixture(random = TRUE, partition = TRUE)
  target <- vecchia_joint_reference(fixture$fit, fixture$newdata, size = 5, method = "covariance")
  set.seed(1)
  draws <- conditional(fixture$fit, fixture$newdata, samples = 4000, output = c("beta", "newdata"),
    local = list(approximation = "vecchia", size = 5, ordering = "none"))
  expect_joint_moments(rbind(draws$beta, draws$newdata), target)
})

test_that("Vecchia respects anisotropy and prediction ordering", {
  fixture <- joint_fixture()
  fit <- fixture$fit
  fit$anisotropy <- TRUE
  fit$coefficients$spcov[["rotate"]] <- 0.7
  fit$coefficients$spcov[["scale"]] <- 0.25
  for (ordering in c("none", "random", "coordinate", "maxmin")) {
    set.seed(1)
    local <- get_local_list_conditional(list(approximation = "vecchia", method = "distance",
      size = 5, ordering = ordering), fit, fixture$newdata)
    coords <- get_conditional_vecchia_covariance(fit, fit$obdata)$coords
    order_o <- conditional_vecchia_order(coords, ordering)
    target <- vecchia_joint_reference(fit, fixture$newdata, size = 5, ord = local$order, order_o = order_o)
    set.seed(1)
    draws <- conditional(fit, fixture$newdata, samples = 4000, output = c("beta", "newdata"),
      local = list(approximation = "vecchia", method = "distance", size = 5, ordering = ordering))
    expect_joint_moments(rbind(draws$beta, draws$newdata), target)
  }
})

test_that("Vecchia local rank handling retains fitted factor contrasts", {
  fixture <- joint_fixture()
  fit <- spglm(y ~ group + offset(off), data = fixture$fit$obdata, family = "poisson",
    xcoord = "cx", ycoord = "cy", contrasts = list(group = "contr.sum"),
    spcov_initial = spcov_initial("exponential", de = 0.4, ie = 0.2, range = 0.3, known = "given"))
  for (method in c("distance", "covariance")) {
    target <- vecchia_joint_reference(fit, fixture$newdata, size = 1, method = method)
    set.seed(1)
    draws <- conditional(fit, fixture$newdata, samples = 4000, output = c("beta", "newdata"),
      local = list(approximation = "vecchia", method = method, size = 1, ordering = "none"))
    expect_joint_moments(rbind(draws$beta, draws$newdata), target)
    expect_identical(rownames(draws$beta), names(coef(fit)))
  }
})

test_that("low-rank conditional moments agree for binomial/small bases", {
  for (family in c("poisson", "binomial")) {
    fixture <- joint_fixture(family = family)
    for (size in c(1, 5)) {
      target <- lowrank_joint_reference(fixture$fit, fixture$newdata, seq_len(size), list(1:2, 3:4))
      set.seed(1)
      draws <- conditional(fixture$fit, fixture$newdata, samples = 5000,
        output = c("beta", "newdata"), local = list(approximation = "low-rank", size_base = size, reorder_base = "none",
          size_new = 2, reorder_new = "none", chunk_size = 300))
      expect_joint_moments(rbind(draws$beta, draws$newdata), target)
    }
  }
})

