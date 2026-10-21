skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set SPMODEL_RUN_EXTRAS = 'true' to run the extras suite")

iid_conditional_fixture <- function(family = NULL, type = "ie", intercept = FALSE, partition = FALSE, local = FALSE, coordinates = FALSE) {
  d <- data.frame(x = seq(-1, 1, length.out = 24), off = seq(-.2, .2, length.out = 24),
    y = rep(c(1, 3, 2, 4, 2, 1), 4), part = rep(c("a", "b"), each = 12))
  if (identical(family, "beta")) d$y <- d$y / 5
  if (identical(family, "binomial")) d$y <- d$y %% 4
  formula <- if (intercept) y ~ 1 + offset(off) else if (identical(family, "binomial"))
    cbind(y, 5 - y) ~ x + offset(off) else y ~ x + offset(off)
  args <- list(formula = formula, data = d, local = local,
    spcov_initial = spcov_initial(type, ie = if (is.null(family)) 2.5 else .5, known = "given"))
  if (partition) args$partition_factor <- ~part
  if (coordinates) {
    args$data$cx <- seq_len(NROW(d))
    args$data$cy <- rep(0:1, 12)
    args$xcoord <- "cx"
    args$ycoord <- "cy"
  }
  if (is.null(family)) {
    args$ddf <- "asymptotic"
    fit <- do.call(splm, args)
  } else {
    args$dispersion_initial <- dispersion_initial(as.character(family),
      dispersion = if (family %in% c("poisson", "binomial")) 1 else 3, known = "given")
    fit <- do.call(spglm, args)
  }
  nd <- data.frame(x = c(-.7, .2, .8), off = c(.1, -.2, .3), part = c("a", "a", "new"))
  if (coordinates) { nd$cx <- 25:27; nd$cy <- c(0, 1, 0) }
  list(fit = fit, newdata = nd)
}

iid_dense_reference <- function(object, newdata, ...) {
  testthat::with_mocked_bindings(
    conditional(object, newdata, ..., local = FALSE),
    get_conditional_iid_variance = function(...) NULL, .package = "spmodel")
}

test_that("Gaussian iid simulation preserves exact draws across chunks and partitions", {
  for (type in c("none", "ie")) for (intercept in c(FALSE, TRUE)) for (partition in c(FALSE, TRUE)) {
    f <- iid_conditional_fixture(type = type, intercept = intercept, partition = partition)
    for (sites in c(1L, 3L)) for (samples in c(1L, 1007L)) {
      nd <- f$newdata[seq_len(sites), , drop = FALSE]
      set.seed(2)
      reference <- iid_dense_reference(f$fit, nd, samples = samples, output = "all")
      set.seed(2)
      actual <- conditional(f$fit, nd, samples = samples, output = "all", local = FALSE)
      expect_identical(actual, reference)
      set.seed(2)
      expect_identical(conditional(f$fit, nd, samples = samples, output = "all"), actual)
    }
  }
  f <- iid_conditional_fixture()
  estimated <- splm(y ~ x + offset(off), f$fit$obdata, spcov_type = "ie", ddf = "asymptotic")
  set.seed(2)
  reference <- iid_dense_reference(estimated, f$newdata, samples = 13)
  set.seed(2)
  expect_identical(conditional(estimated, f$newdata, samples = 13), reference)
  f <- iid_conditional_fixture(local = list(index = rep(1:2, each = 12), var_adjust = "theoretical"))
  set.seed(2)
  reference <- iid_dense_reference(f$fit, f$newdata, samples = 13)
  set.seed(2)
  expect_identical(conditional(f$fit, f$newdata, samples = 13), reference)
})

test_that("iid GLM draws retain the dense joint moments and coefficient stream", {
  for (family in c("poisson", "binomial", "nbinomial", "Gamma", "inverse.gaussian", "beta")) {
    for (type in c("none", "ie")) {
      f <- iid_conditional_fixture(family, type, partition = TRUE)
      fit <- f$fit
      n <- fit$n
      X <- model.matrix(fit)
      X0 <- model.matrix(~x, f$newdata)
      B <- vcov(fit)
      v <- if (type == "none") 1e-4 else .5
      eta <- fitted(fit, type = "link")
      D <- as.matrix(get_D(family, eta, fit$y, fit$size, as.vector(coef(fit, "dispersion"))))
      Q <- solve(diag(1 / v, n) - D)
      A <- Q %*% (X / v)
      covariance <- rbind(cbind(B, B %*% t(A), B %*% t(X0)),
        cbind(A %*% B, Q + A %*% B %*% t(A), A %*% B %*% t(X0)),
        cbind(X0 %*% B, X0 %*% B %*% t(A), v * diag(3) + X0 %*% B %*% t(X0)))
      center <- get_latent_center(fit)
      mu <- c(center$betahat, center$eta, X0 %*% center$betahat + f$newdata$off)
      set.seed(2)
      reference <- iid_dense_reference(fit, f$newdata, samples = 1007, output = c("beta", "latent", "newdata"))
      set.seed(2)
      actual <- conditional(fit, f$newdata, samples = 1007, output = c("beta", "latent", "newdata"))
      expect_identical(actual$beta, reference$beta)
      expect_identical(rownames(actual$latent), rownames(X))
      expect_identical(rownames(actual$newdata), rownames(X0))
      set.seed(2)
      expect_identical(conditional(fit, f$newdata, samples = 1007,
        output = c("beta", "latent", "newdata"), local = FALSE), actual)
      set.seed(2)
      draws <- conditional(fit, f$newdata, samples = 16000, output = c("beta", "latent", "newdata"))
      values <- do.call(rbind, draws)
      se <- sqrt((outer(diag(covariance), diag(covariance)) + covariance^2) / (NCOL(values) - 1))
      expect_true(all(abs(rowMeans(values) - mu) < 7 * sqrt(diag(covariance) / NCOL(values))))
      expect_true(all(abs(cov(t(values)) - covariance) < 7 * se))
      set.seed(2)
      expect_identical(conditional(fit, output = "beta", samples = 1007), actual$beta)
      expect_identical(dim(conditional(fit, output = "latent", samples = 1)), c(n, 1L))
      expect_equal(conditional(fit, output = "object", samples = 1)[, 1], unname(eta))
    }
  }
})

test_that("iid simulation skips dense preparation and unused latent draws", {
  for (family in list(NULL, "poisson")) {
    f <- iid_conditional_fixture(family, partition = TRUE)
    nd <- f$newdata[rep(1:3, length.out = 5001), ]
    fail <- function(...) stop("unexpected dense or spatial preparation")
    testthat::with_mocked_bindings({
      expect_message(actual <- conditional(f$fit, nd, samples = 3), NA)
      expect_identical(dim(actual), c(5001L, 3L))
      expect_identical(dim(conditional(f$fit, output = "beta", samples = 1)), c(2L, 1L))
    }, covmatrix = fail, spdist = fail, spdist_vectors = fail,
      get_conditional_glm_joint = fail, get_conditional_glm_base = fail,
      get_conditional_new_from_base_adjust = fail, prepare_conditional_vecchia = fail,
      get_D = fail, get_extra_partition_list = fail, model_matrix_group_labels = fail,
      .package = "spmodel")
  }
  d <- data.frame(x = seq(-1, 1, length.out = 5001), y = rep(c(1, 3, 2), 1667))
  fit <- splm(y ~ x, d, spcov_type = "ie", ddf = "asymptotic")
  testthat::with_mocked_bindings({
    expect_message(value <- conditional(fit, d[1:3, ], samples = 2), NA)
    expect_identical(dim(value), c(3L, 2L))
  }, covmatrix = function(...) stop("unexpected covariance"), .package = "spmodel")
})

test_that("iid eligibility preserves model exclusions and zero-variance behavior", {
  f <- iid_conditional_fixture()
  for (n in c(4999, 5000, 5001)) for (m in c(4999, 5000, 5001)) {
    fit <- f$fit
    fit$n <- n
    testthat::with_mocked_bindings({
      expect_message(conditional(fit, f$newdata[rep(1, m), ], samples = 1), NA)
    }, draw_conditional_iid_error = function(n, samples, variance, ...) matrix(0, n, samples), .package = "spmodel")
  }
  expect_null(get_conditional_iid_variance(f$fit, simulate_covparams = TRUE))
  fit <- f$fit
  fit$random <- ~g
  expect_null(get_conditional_iid_variance(fit))
  fit <- f$fit
  class(fit) <- "spautor"
  expect_null(get_conditional_iid_variance(fit))
  fit <- f$fit
  fit$coefficients$spcov <- spcov_params("exponential", de = 0, ie = 2.5, range = 1)
  expect_null(get_conditional_iid_variance(fit))
  fit <- f$fit
  fit$coefficients$spcov[["ie"]] <- 0
  expect_null(get_conditional_iid_variance(fit))
  expect_error(conditional(fit, f$newdata, samples = 1, local = FALSE), "positive")
  for (v in c(Inf, NA_real_)) {
    fit$coefficients$spcov[["ie"]] <- v
    expect_error(conditional(fit, f$newdata, samples = 1, local = FALSE), "variance")
  }
  testthat::with_mocked_bindings({
    expect_error(get_conditional_iid_variance(f$fit), "variance")
  }, spcov_ie_stabilized = function(...) -1, .package = "spmodel")
  expect_error(conditional(f$fit, f$newdata, output = "random"), "random effects")
  for (family in list(NULL, "poisson")) {
    f <- iid_conditional_fixture(family, partition = TRUE)
    nd <- f$newdata[c("x", "off")]
    expect_error(conditional(f$fit, nd, samples = 1), "part.*not found")
    expect_error(iid_dense_reference(f$fit, nd, samples = 1), "part.*not found")
  }
})

test_that("iid latent precision checks cover positive beta curvature", {
  f <- iid_conditional_fixture("beta")
  for (d in c(.1, .5, 1, Inf, NA_real_)) {
    testthat::with_mocked_bindings({
      if (identical(d, .1)) {
        expect_true(all(is.finite(get_conditional_iid_glm(f$fit, matrix(numeric(), 0, 2),
          variance = 2, samples = 3, latent = TRUE)$latent)))
      } else {
        expect_error(get_conditional_iid_glm(f$fit, matrix(numeric(), 0, 2),
          variance = 2, samples = 3, latent = TRUE), "latent covariance is not positive definite")
      }
    }, get_D = function(...) Matrix::Diagonal(f$fit$n, d), .package = "spmodel")
  }
})

iid_ignored_local <- list(NULL, FALSE, TRUE, list(approximation = "vecchia", size = 3),
  list(approximation = "low-rank", size_base = 8, size_new = 2, reorder_base = "grts", kmeans_new = TRUE),
  list(approximation = "bad", size_random = 0, chunk_size = 0, parallel = TRUE, ncores = -1), NA, "unused")

test_that("iid conditional simulation ignores every local setting before spatial preparation", {
  for (family in list(NULL, "poisson")) for (type in c("none", "ie")) for (coordinates in c(FALSE, TRUE)) {
    f <- iid_conditional_fixture(family, type, coordinates = coordinates)
    output <- c("newdata", "beta", if (!is.null(family)) "latent")
    set.seed(18)
    reference <- conditional(f$fit, f$newdata, samples = 7, output = output, local = FALSE)
    fail <- function(...) stop("unexpected spatial preparation")
    resolve_lowrank <- get_local_list_conditional_lowrank
    testthat::with_mocked_bindings({
      for (ctl in iid_ignored_local) {
        set.seed(18)
        expect_message(actual <- conditional(f$fit, f$newdata, samples = 7, output = output, local = ctl), NA)
        expect_identical(actual, reference)
      }
    }, covmatrix = fail, spdist = fail, get_local_list_conditional_lowrank = function(local, ...) {
      expect_identical(local$method_base, "all")
      expect_identical(local$method_new, "all")
      resolve_lowrank(local, ...)
    },
      get_local_list_conditional_vecchia = fail, prepare_conditional_vecchia = fail,
      .package = "spmodel")
    if (is.null(family)) {
      reference <- loocv(f$fit, cv_predict = TRUE, se.fit = TRUE, interval = "prediction", local = FALSE)
      testthat::with_mocked_bindings({
        for (ctl in iid_ignored_local) {
          expect_identical(loocv(f$fit, cv_predict = TRUE, se.fit = TRUE,
            interval = "prediction", local = ctl), reference)
        }
      }, get_local_list_prediction = fail, covmatrix = fail, .package = "spmodel")
    }
  }
})

test_that("iid local override retains exact covariance-parameter simulation and zero errors", {
  f <- iid_conditional_fixture()
  fit <- splm(y ~ x + offset(off), f$fit$obdata, spcov_type = "ie", ddf = "satterthwaite")
  set.seed(18)
  reference <- conditional(fit, f$newdata, samples = 7, output = c("newdata", "beta", "spcov"),
    simulate_covparams = TRUE, local = FALSE)
  expect_gt(stats::var(reference$spcov["ie", ]), 0)
  for (ctl in iid_ignored_local) {
    set.seed(18)
    expect_message(actual <- conditional(fit, f$newdata, samples = 7,
      output = c("newdata", "beta", "spcov"), simulate_covparams = TRUE, local = ctl), NA)
    expect_identical(actual, reference)
  }
  fit$coefficients$spcov[["ie"]] <- 0
  for (ctl in iid_ignored_local) {
    expect_error(conditional(fit, f$newdata, samples = 1, local = ctl), "positive")
  }
})

test_that("models with spatial covariance or random effects retain local controls", {
  for (family in list(NULL, "poisson")) for (random in c(FALSE, TRUE)) {
    f <- iid_conditional_fixture(family, coordinates = TRUE)
    if (random) f$fit$random <- ~g else
      f$fit$coefficients$spcov <- spcov_params("exponential", de = 1, ie = .5, range = 1)
    ctl <- list(approximation = "bad")
    testthat::with_mocked_bindings({
      expect_error(conditional(f$fit, f$newdata, samples = 1, local = ctl), "retained local controls")
    }, get_local_list_conditional = function(local, ...) {
      expect_identical(local, ctl)
      stop("retained local controls")
    }, .package = "spmodel")
  }
})
