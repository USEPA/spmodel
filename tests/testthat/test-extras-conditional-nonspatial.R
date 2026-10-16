skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set SPMODEL_RUN_EXTRAS = 'true' to run the extras suite")

nonspatial_fixture <- function(family = NULL, type = "none", coordinates = FALSE,
                               random = FALSE, partition = FALSE, estmethod = "reml") {
  d <- data.frame(x = seq(-1, 1, length.out = 24), off = seq(-0.2, 0.2, length.out = 24),
    cx = seq_len(24), cy = rep(c(0, 1), 12), g = rep(letters[1:4], 6),
    part = rep(c("a", "b"), each = 12), y = rep(c(1, 3, 2, 4, 2, 1), 4))
  nd <- data.frame(x = c(-0.7, 0.2, 0.8), off = c(0.1, -0.2, 0.3),
    cx = c(30, 30, 32), cy = c(0, 0, 1), g = c("a", "a", "new"), part = c("a", "a", "new"))
  if (identical(family, "beta")) d$y <- d$y / 5
  if (identical(family, "binomial")) d$y <- d$y %% 4
  if (!coordinates) {
    d[c("cx", "cy")] <- NULL
    nd[c("cx", "cy")] <- NULL
  }
  args <- list(formula = if (identical(family, "binomial")) cbind(y, 5 - y) ~ x + offset(off) else y ~ x + offset(off),
    data = d, local = FALSE, estmethod = estmethod,
    spcov_initial = spcov_initial(type, ie = if (is.null(family)) 2.5 else 0.5,
      de = if (type == "exponential") 0 else NULL,
      range = if (type == "exponential") 1 else NULL, known = "given"))
  if (coordinates) { args$xcoord <- "cx"; args$ycoord <- "cy" }
  if (random) {
    args$random <- ~ (1 | g) + (0 + x | g)
    args$randcov_initial <- randcov_initial(g = 0.3, "0 + x | g" = 0.1, known = "given")
  }
  if (partition) args$partition_factor <- ~part
  if (is.null(family)) {
    args$ddf <- "asymptotic"
    fit <- do.call(splm, args)
  } else {
    args$dispersion_initial <- dispersion_initial(as.character(family),
      dispersion = if (family %in% c("poisson", "binomial")) 1 else 3, known = "given")
    fit <- do.call(spglm, args)
  }
  list(fit = fit, data = d, newdata = nd)
}

expect_normal_moments <- function(draws, mean, covariance) {
  n <- NCOL(draws)
  expect_true(all(abs(rowMeans(draws) - mean) < 6 * sqrt(diag(covariance) / n)))
  se <- sqrt((outer(diag(covariance), diag(covariance)) + covariance^2) / (n - 1))
  expect_true(all(abs(cov(t(draws)) - covariance) < 6 * se))
}

test_that("independent Gaussian draws retain nugget and coefficient uncertainty", {
  for (type in c("exponential", "none", "ie")) {
    f <- nonspatial_fixture(type = type, coordinates = TRUE)
    expect_equal(unname(coef(f$fit, "spcov")[["ie"]]), 2.5)
    for (base in c(8L, 24L)) {
      settings <- list(approximation = "low-rank", size_base = base, size_new = 2,
        reorder_base = "none", reorder_new = "none", chunk_size = 7)
      for (ctl in list(FALSE, settings)) {
        for (sites in c(1L, 3L)) {
          nd <- f$newdata[seq_len(sites), , drop = FALSE]
          for (samples in c(1L, 15L)) {
            set.seed(2)
            expect_identical(dim(conditional(f$fit, nd, samples = samples, local = ctl)), c(sites, samples))
          }
        }
        X <- model.matrix(~x, f$newdata)
        B <- vcov(f$fit)
        set.seed(2)
        draws <- conditional(f$fit, f$newdata, samples = 12000, local = ctl,
          output = c("newdata", "beta"))
        expected <- 2.5 * diag(3) + X %*% B %*% t(X)
        expect_normal_moments(draws$newdata, as.numeric(X %*% coef(f$fit)) + f$newdata$off, expected)
        residual <- draws$newdata - X %*% draws$beta - f$newdata$off
        expect_normal_moments(residual, rep(0, 3), 2.5 * diag(3))
      }
    }
    set.seed(2)
    z <- matrix(rnorm(9), 3, 3)
    set.seed(2)
    value <- get_conditional_new_from_base_adjust(list(newdata = f$newdata), f$fit,
      matrix(0, 24, 1), matrix(0, 24, 2), matrix(0, 2, 3), sqrt(2.5) * diag(24), 3)
    L <- value %*% solve(z)
    expect_equal(tcrossprod(L), 2.5 * diag(3), tolerance = 1e-12)
    expect_true(is.matrix(value))
  }
})

test_that("nonspatial covariance and simulation do not need geometric distances", {
  for (family in list(NULL, "poisson")) for (type in c("none", "ie")) {
    f <- nonspatial_fixture(family, type)
    with_coords <- nonspatial_fixture(family, type, coordinates = TRUE)
    expect_equal(coef(f$fit), coef(with_coords$fit))
    expect_equal(predict(f$fit, f$newdata, se.fit = TRUE, local = FALSE),
      predict(with_coords$fit, with_coords$newdata, se.fit = TRUE, local = FALSE))
    v <- if (is.null(family)) 2.5 else if (type == "none") 1e-4 else 0.5
    expect_equal(unname(coef(f$fit, "spcov")[["ie"]]), if (!is.null(family) && type == "none") 0 else v)
    testthat::with_mocked_bindings({
      for (sites in c(1L, 3L)) {
        nd <- f$newdata[seq_len(sites), , drop = FALSE]
        expect_equal(covmatrix(f$fit), v * diag(24))
        expect_equal(covmatrix(f$fit, nd, cov_type = "pred.pred"), v * diag(sites))
        expect_equal(covmatrix(f$fit, nd, cov_type = "pred.obs"), matrix(0, sites, 24))
        expect_equal(covmatrix(f$fit, nd, cov_type = "obs.pred"), matrix(0, 24, sites))
        expect_length(predict(f$fit, nd, se.fit = TRUE, local = FALSE)$fit, sites)
        for (samples in c(1L, 11L)) {
          set.seed(2)
          expect_identical(dim(conditional(f$fit, nd, samples = samples, local = FALSE)), c(sites, samples))
        }
      }
    }, spdist = function(...) stop("unexpected distance calculation"),
       spdist_vectors = function(...) stop("unexpected cross-distance calculation"), .package = "spmodel")
    set.seed(2)
    a <- conditional(f$fit, f$newdata, samples = 11, local = FALSE)
    set.seed(2)
    expect_equal(a, conditional(with_coords$fit, with_coords$newdata, samples = 11, local = FALSE))
    set.seed(2)
    expect_identical(dim(conditional(f$fit, output = "beta", samples = 1, local = FALSE)), c(2L, 1L))
    for (base in c(8L, 24L)) {
      set.seed(2)
      expect_identical(dim(conditional(f$fit, f$newdata, samples = 11,
        local = list(approximation = "low-rank", size_base = base, size_new = 2,
          reorder_base = "none", reorder_new = "none", chunk_size = 4))), c(3L, 11L))
    }
    expect_identical(dim(conditional(f$fit, f$newdata, samples = 1, local = TRUE)), c(3L, 1L))
    expect_identical(dim(conditional(f$fit, f$newdata, samples = 1,
      local = list(approximation = "low-rank", size_base = 8))), c(3L, 1L))
    expect_equal(dim(conditional(f$fit, f$newdata, samples = 1, local = list(method = "all"))), c(3L, 1L))
  }
})

test_that("nonspatial families preserve link transformations and response support", {
  for (family in c("poisson", "binomial", "nbinomial", "Gamma", "inverse.gaussian", "beta")) {
    for (type in c("none", "ie")) {
      f <- nonspatial_fixture(family, type)
      for (sites in c(1L, 3L)) for (samples in c(1L, 13L)) {
        nd <- f$newdata[seq_len(sites), , drop = FALSE]
        sizes <- seq_len(sites) + 3L
        set.seed(2)
        link <- conditional(f$fit, nd, samples = samples, local = FALSE, type = "link", newdata_size = sizes)
        set.seed(2)
        response <- conditional(f$fit, nd, samples = samples, local = FALSE, type = "response", newdata_size = sizes)
        expected <- if (family %in% c("binomial", "beta")) plogis(link) else exp(link)
        if (family == "binomial") expected <- expected * sizes
        expect_equal(response, expected)
        set.seed(2)
        fresh <- conditional(f$fit, nd, samples = samples, local = FALSE, type = "new", newdata_size = sizes)
        expect_identical(dim(fresh), c(sites, samples))
        expect_true(all(is.finite(fresh) & fresh >= 0))
        if (family %in% c("poisson", "binomial", "nbinomial")) expect_equal(fresh, round(fresh))
        if (family == "binomial") expect_true(all(fresh <= sizes))
        if (family == "beta") expect_true(all(fresh > 0 & fresh < 1))
      }
    }
  }
})

test_that("nonspatial mixed-model covariance preserves groups and partitions", {
  for (family in list(NULL, "poisson")) for (type in c("none", "ie")) {
    f <- nonspatial_fixture(family, type, random = TRUE, partition = TRUE)
    d <- f$data
    nd <- f$newdata
    v <- if (is.null(family)) 2.5 else if (type == "none") 1e-4 else 0.5
    rcov <- function(a, b) (0.3 + 0.1 * outer(a$x, b$x)) *
      outer(a$g, b$g, "==") * outer(a$part, b$part, "==")
    V <- v * diag(24) + rcov(d, d)
    C <- rcov(nd, d)
    T <- v * diag(3) + rcov(nd, nd)
    expect_equal(unname(covmatrix(f$fit)), V)
    expect_equal(unname(covmatrix(f$fit, nd, cov_type = "pred.pred")), T)
    expect_equal(unname(covmatrix(f$fit, nd, cov_type = "pred.obs")), C)
    expect_equal(unname(covmatrix(f$fit, nd, cov_type = "obs.pred")), t(C))
    expect_gt(T[1, 2], 0)
    split <- nd
    split$part[2] <- "b"
    expect_equal(covmatrix(f$fit, split, cov_type = "pred.pred")[1, 2], 0)
    X <- model.matrix(~x, d)
    Xn <- model.matrix(~x, nd)
    B <- vcov(f$fit)
    K <- C %*% solve(V)
    beta <- coef(f$fit)
    if (is.null(family)) {
      H <- Xn - K %*% X
      mu <- as.numeric(Xn %*% beta + K %*% (d$y - d$off - X %*% beta)) + nd$off
      covariance <- T - K %*% t(C) + H %*% B %*% t(H)
    } else {
      eta <- fitted(f$fit, type = "link")
      Q <- solve(solve(V) + diag(exp(eta)))
      A <- Q %*% solve(V, X)
      H <- Xn + K %*% (A - X)
      mu <- as.numeric(Xn %*% beta + K %*% (eta - d$off - X %*% beta)) + nd$off
      covariance <- T - K %*% t(C) + K %*% Q %*% t(K) + H %*% B %*% t(H)
    }
    set.seed(2)
    draws <- conditional(f$fit, nd, samples = 16000, local = FALSE, output = c("newdata", "random"))
    expect_normal_moments(draws$newdata, mu, covariance)
    expect_true(any(!attr(draws$random, "effect_info")$conditioned))
  }
})

test_that("independent latent draws preserve the effective variance floor", {
  for (type in c("none", "ie")) {
    f <- nonspatial_fixture("poisson", type)
    X <- model.matrix(~x, f$newdata)
    v <- if (type == "none") 1e-4 else 0.5
    set.seed(2)
    draws <- conditional(f$fit, f$newdata, samples = 16000, local = FALSE, output = c("newdata", "beta"))
    residual <- draws$newdata - X %*% draws$beta - f$newdata$off
    expect_normal_moments(residual, rep(0, 3), v * diag(3))
    expect_normal_moments(draws$newdata, as.numeric(X %*% coef(f$fit)) + f$newdata$off,
      X %*% vcov(f$fit) %*% t(X) + v * diag(3))
  }
})

test_that("nonspatial simulation accepts geometry and implicit prediction data", {
  for (family in list(NULL, "poisson")) {
    f <- nonspatial_fixture(family, "ie", coordinates = TRUE)
    sfnew <- sf::st_as_sf(f$newdata, coords = c("cx", "cy"))
    set.seed(2)
    expected <- conditional(f$fit, f$newdata, samples = 5, local = FALSE)
    set.seed(2)
    expect_equal(conditional(f$fit, sfnew, samples = 5, local = FALSE), expected)
    f$fit$newdata <- f$newdata[c("x", "off", "g", "part")]
    set.seed(2)
    expect_equal(conditional(f$fit, samples = 5, local = FALSE), expected)
  }
})

test_that("coordinate-free fitting has an explicit estimation-method boundary", {
  for (type in c("none", "ie")) {
    for (method in c("reml", "ml")) expect_s3_class(nonspatial_fixture(type = type, estmethod = method)$fit, "splm")
    for (method in c("sv-wls", "sv-cl")) expect_error(nonspatial_fixture(type = type, estmethod = method), "semivariogram estimation requires coordinates")
  }
  f <- nonspatial_fixture(type = "exponential", coordinates = TRUE)
  expect_error(conditional(f$fit, f$newdata[c("x", "off")], local = FALSE), "Coordinate column")
})

test_that("large random-effect designs warn and continue with/without observation blocking", {
  testthat::local_mocked_bindings(conditional_spatial_object = function(...) stop("past size warning"), .package = "spmodel")
  for (q in c(10001L, 20000L)) {
    design <- list(observed = Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(1L, q)),
      variance = rep(1, q), info = data.frame(conditioned = rep(TRUE, q)))
    for (ctl in list(NULL, list(approximation = "vecchia"))) {
      expect_warning(
        expect_error(prepare_conditional_random(list(), design, ctl), "past size warning"),
        paste0(format(q, big.mark = ",", scientific = FALSE), " random-effect coefficients.*exceedingly long computational times"))
    }
  }
  design <- list(observed = Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(1L, 10001L)),
    variance = rep(1, 10001), info = data.frame(conditioned = rep(TRUE, 10001)))
  design$variance[10001] <- 0
  expect_warning(expect_error(prepare_conditional_random(list(), design), "past size warning"), NA)
  design$variance[10001] <- 1
  design$info$conditioned[10001] <- FALSE
  expect_warning(expect_error(prepare_conditional_random(list(), design), "past size warning"), NA)
})
