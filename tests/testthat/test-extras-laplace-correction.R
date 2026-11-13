skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

cv_prediction_fixture <- function(family = "poisson", areal = FALSE, partition = FALSE, random = FALSE) {
  set.seed(203)
  n <- 24
  d <- data.frame(x = rnorm(n), cx = runif(n), cy = runif(n), off = runif(n, -.2, .2),
    group = rep(letters[1:4], 6), part = rep(c("a", "b"), each = 12), trials = rep(4:7, 6))
  d$y <- switch(family, poisson = rpois(n, exp(.5 + .2 * d$x + d$off)),
    nbinomial = rnbinom(n, mu = exp(.5 + .2 * d$x + d$off), size = 2),
    binomial = rbinom(n, d$trials, plogis(.2 * d$x + d$off)))
  args <- list(formula = if (family == "binomial") cbind(y, trials - y) ~ x + offset(off) else y ~ x + offset(off),
    data = d, dispersion_initial = dispersion_initial(as.character(family), dispersion = if (family == "nbinomial") 2 else 1, known = "given"))
  if (partition) args$partition_factor <- ~part
  if (random) {
    args$random <- ~group
    args$randcov_initial <- randcov_initial(group = .15, known = "given")
  }
  if (areal) {
    args$data$y[22:24] <- NA
    args$W <- 1 * (abs(outer(1:n, 1:n, "-")) == 1)
    args$spcov_initial <- spcov_initial("car", de = .3, ie = .1, range = .2, known = "given")
    fit <- do.call(spgautor, args)
  } else {
    args$xcoord <- "cx"; args$ycoord <- "cy"
    args$spcov_initial <- spcov_initial("exponential", de = .3, ie = .1, range = .2, known = "given")
    fit <- do.call(spglm, args)
  }
  list(fit = fit, args = args)
}

test_that("nonlocal GLM CV agrees with held-out kriging weights", {
  for (family in c("poisson", "nbinomial", "binomial")) for (areal in c(FALSE, TRUE)) {
    fit <- cv_prediction_fixture(family, areal)$fit
    X <- model.matrix(fit); S <- as.matrix(covmatrix(fit)); eta <- fitted(fit, "link")
    offset <- model.offset(model.frame(fit)); w <- eta - offset
    for (folds in list(seq_len(fit$n), rep(1:3, length.out = fit$n))) {
      expected <- numeric(fit$n)
      for (F in split(seq_len(fit$n), folds)) {
        Q <- solve(S[-F, -F, drop = FALSE]); Xm <- X[-F, , drop = FALSE]
        B <- solve(crossprod(Xm, Q %*% Xm), crossprod(Xm, Q))
        L <- X[F, , drop = FALSE] %*% B + S[F, -F, drop = FALSE] %*% Q %*% (diag(nrow(Xm)) - Xm %*% B)
        expected[F] <- as.numeric(L %*% (w[-F])) + offset[F]
      }
      args <- list(object = fit, folds_index = folds, cv_predict = TRUE, se.fit = TRUE, local = FALSE)
      on <- do.call(kcv, args)
      expect_equal(on$cv_predict, expected, tolerance = 1e-9)
      response <- invlink(expected, family, fit$size)
      expect_equal(on$stats$bias, mean(fit$y - response))
      expect_equal(on$stats$MSPE, mean((fit$y - response)^2))
      expect_equal(on$stats$RMSPE, sqrt(mean((fit$y - response)^2)))
      rp <- kcv(fit, folds_index = folds, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE)
      derivative <- if (family == "binomial") fit$size * plogis(expected) * (1 - plogis(expected)) else exp(expected)
      expect_equal(rp$cv_predict, response)
      expect_equal(rp$se.fit, derivative * on$se.fit)
      if (length(unique(folds)) == fit$n) {
        loo <- loocv(fit, cv_predict = TRUE, se.fit = TRUE)
        expect_equal(loo, on, tolerance = 1e-9)
        expect_identical(kcv(fit, k = fit$n, cv_predict = TRUE, se.fit = TRUE), loo)
      }
    }

  }
})

test_that("local GLM LOO retains fitted centers and mode curvature", {
  fit <- cv_prediction_fixture("binomial", partition = TRUE, random = TRUE)$fit
  X <- model.matrix(fit); S <- as.matrix(covmatrix(fit)); d <- fit$obdata
  w <- fitted(fit, "link") - d$off
  beta <- coef(fit)
  for (method in c("distance", "covariance")) {
    ctl <- list(method = method, size = 8)
    expected <- vapply(seq_len(fit$n), function(i) {
      pool <- which(seq_len(fit$n) != i & d$part == d$part[i])
      score <- if (method == "distance") (d$cx[pool] - d$cx[i])^2 + (d$cy[pool] - d$cy[i])^2 else -abs(S[i, pool])
      N <- pool[order(score)[seq_len(min(8, length(pool)))]]
      as.numeric(X[i, , drop = FALSE] %*% beta + S[i, N, drop = FALSE] %*%
        solve(S[N, N, drop = FALSE], w[N] - X[N, , drop = FALSE] %*% beta)) + d$off[i]
    }, numeric(1))
    on <- loocv(fit, local = ctl, cv_predict = TRUE, se.fit = TRUE)
    expect_equal(on$cv_predict, expected, tolerance = 1e-9)
    response <- loocv(fit, local = ctl, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE)
    prob <- plogis(on$cv_predict)
    expect_equal(response$cv_predict, invlink(on$cv_predict, fit$family, fit$size))
    expect_equal(response$se.fit, fit$size * prob * (1 - prob) * on$se.fit)
  }
})

test_that("local kcv matches training-fold predictions with full-fit coefficient covariance", {
  z <- cv_prediction_fixture("binomial", partition = TRUE, random = TRUE)
  z$args$local <- list(index = rep(1:3, 8), var_adjust = "none")
  fit <- do.call(spglm, z$args)
  folds <- rep(1:3, each = 8)
  ctl <- list(method = "distance", size = 8)
  {
    expected <- se <- numeric(fit$n)
    for (F in split(seq_len(fit$n), folds)) {
      args <- z$args
      args$data$y[F] <- NA
      args$local$index <- fit$local_index[-F]
      train <- do.call(spglm, args)
      train$vcov$fixed <- fit$vcov$fixed
      pred <- predict(train, se.fit = TRUE, local = ctl)
      expected[F] <- pred$fit; se[F] <- pred$se.fit
    }
    actual <- kcv(fit, folds_index = folds, local = ctl, cv_predict = TRUE, se.fit = TRUE)
    expect_equal(actual$cv_predict, expected, tolerance = 1e-9)
    expect_equal(actual$se.fit, se, tolerance = 1e-9)
  }
})

test_that("GLM CV agrees between serial and parallel paths", {
  for (areal in c(FALSE, TRUE)) {
    z <- cv_prediction_fixture(areal = areal)
    fit <- z$fit
    if (!areal) {
      z$args$local <- list(index = rep(1:3, 8), var_adjust = "none")
      fit <- do.call(spglm, z$args)
    }
    methods <- if (areal) "all" else c("all", "distance")
    for (method in methods) {
      args <- list(object = fit, cv_predict = TRUE, se.fit = TRUE,
        local = list(method = method, size = 8, parallel = FALSE))
      parallel_args <- args
      parallel_args$local$parallel <- TRUE
      parallel_args$local$ncores <- 2
      expect_equal(do.call(loocv, parallel_args), do.call(loocv, args), tolerance = 1e-12)
      args$folds_index <- parallel_args$folds_index <- rep(1:3, length.out = fit$n)
      expect_equal(do.call(kcv, parallel_args), do.call(kcv, args), tolerance = 1e-12)
    }
  }
})
