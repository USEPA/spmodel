skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

laplace_fixture <- function(family = "poisson", estmethod = "reml", local = FALSE, areal = FALSE, iid = FALSE) {
  set.seed(1910)
  n <- 30
  d <- data.frame(x = rnorm(n), east = runif(n), north = runif(n), off = runif(n, -.2, .2))
  mu <- exp(1 + .2 * d$x + d$off)
  d$y <- switch(family, poisson = rpois(n, mu), nbinomial = rnbinom(n, mu = mu, size = 2),
    binomial = rbinom(n, 8, plogis(.2 * d$x + d$off)))
  formula <- if (family == "binomial") cbind(y, 8 - y) ~ x + offset(off) else y ~ x + offset(off)
  args <- list(formula = formula, data = d, family = family, estmethod = estmethod)
  if (family == "nbinomial") {
    args$family <- NULL
    args$dispersion_initial <- dispersion_initial("nbinomial", dispersion = 2, known = "given")
  }
  if (areal) {
    args$data$y[27:30] <- NA
    args$W <- 1 * (abs(outer(seq_len(n), seq_len(n), "-")) == 1)
    args$spcov_initial <- spcov_initial("car", de = .3, ie = .1, range = .2, known = "given")
    fit <- do.call(spgautor, args)
  } else {
    args$spcov_initial <- if (iid) spcov_initial("ie", ie = .3, known = "given") else
      spcov_initial("exponential", de = .3, ie = .1, range = .2, known = "given")
    args$xcoord <- "east"; args$ycoord <- "north"; args$local <- local
    fit <- do.call(spglm, args)
  }
  new <- d[1:4, ]; new$east <- new$east + .03
  list(fit = fit, new = new)
}

test_that("third derivatives agree with central differences for every family", {
  for (family in c("poisson", "nbinomial", "binomial", "Gamma", "inverse.gaussian", "beta")) {
    grid <- expand.grid(w = c(-4, -1, 0, .7, 3), y = if (family == "beta") c(.02, .3, .8, .98) else c(.2, 1, 5),
      size = c(1, 8), dispersion = c(.3, 2, 15))
    numeric_D3 <- (diag(get_D(family, grid$w + 1e-5, grid$y, grid$size, grid$dispersion)) -
      diag(get_D(family, grid$w - 1e-5, grid$y, grid$size, grid$dispersion))) / 2e-5
    actual <- get_D3(family, grid$w, grid$y, grid$size, grid$dispersion)
    expect_type(actual, "double")
    expect_lt(max(abs(actual - numeric_D3) / pmax(1, abs(numeric_D3))), 1e-5)
  }
})

test_that("stored corrections use projected curvature and original observation order", {
  for (family in c("poisson", "nbinomial", "binomial")) {
    for (method in c("ml", "reml")) {
      fit <- laplace_fixture(family, method)$fit
      X <- model.matrix(fit)
      Si <- solve(as.matrix(covmatrix(fit)))
      B <- solve(crossprod(X, Si %*% X), crossprod(X, Si))
      P <- Si - Si %*% X %*% B
      eta <- fitted(fit, type = "link")
      phi <- as.vector(coef(fit, type = "dispersion"))
      D <- get_D(family, eta, fit$y, fit$size, phi)
      V <- solve(P - as.matrix(D))
      t3 <- get_D3(family, eta, fit$y, fit$size, phi)
      delta <- .5 * as.numeric(V %*% (t3 * diag(V)))
      expect_equal(unname(fit$laplace$delta), delta, tolerance = 1e-4)
      expect_identical(names(fit$laplace$delta), names(eta))
      expect_equal(unname(fit$laplace$betahat), as.numeric(B %*%
        (eta - model.offset(model.frame(fit)) + fit$laplace$delta)), tolerance = 1e-10)
      if (method == "ml") {
        fixed <- solve(Si - as.matrix(D))
        expect_gt(max(abs(delta - .5 * as.numeric(fixed %*% (t3 * diag(fixed))))), 1e-4)
      }
    }
  }
  for (areal in c(FALSE, TRUE)) {
    fit <- laplace_fixture(local = if (areal) FALSE else list(method = "random", size = 10), areal = areal)$fit
    Sigma <- as.matrix(covmatrix(fit))
    if (!areal) {
      expect_false(identical(order(fit$local_index), seq_len(fit$n)))
      Sigma[outer(fit$local_index, fit$local_index, "!=")] <- 0
    }
    X <- model.matrix(fit); Si <- solve(Sigma)
    B <- solve(crossprod(X, Si %*% X), crossprod(X, Si))
    eta <- fitted(fit, type = "link")
    V <- solve(Si - Si %*% X %*% B + diag(exp(eta)))
    expected <- -.5 * as.numeric(V %*% (exp(eta) * diag(V)))
    expect_equal(unname(fit$laplace$delta), expected, tolerance = 1e-4)
    expect_equal(unname(fit$laplace$betahat), as.numeric(B %*%
      (eta - model.offset(model.frame(fit)) + fit$laplace$delta)), tolerance = 1e-10)
  }
})

test_that("prediction shifts only the center and old objects retain mode predictions", {
  for (areal in c(FALSE, TRUE)) {
    fixture <- laplace_fixture(areal = areal)
    fit <- fixture$fit
    args <- if (areal) list(object = fit) else list(object = fit, newdata = fixture$new)
    old <- fit; old$laplace <- NULL
    for (interval in c("none", "prediction", "confidence")) {
      off <- do.call(predict, c(args, list(se.fit = TRUE, interval = interval, fit_correct = FALSE)))
      on <- do.call(predict, c(args, list(se.fit = TRUE, interval = interval)))
      oldargs <- args; oldargs$object <- old
      expect_identical(do.call(predict, c(oldargs, list(se.fit = TRUE, interval = interval))), off)
      expect_identical(on$se.fit, off$se.fit)
      expect_gt(max(abs(on$fit - off$fit)), 0)
    }
    weights <- do.call(predict, c(args, list(type = "weight")))
    expect_identical(weights, do.call(predict, c(args, list(type = "weight", fit_correct = FALSE))))
    expect_identical(do.call(predict, c(args, list(type = "terms"))),
      do.call(predict, c(args, list(type = "terms", fit_correct = FALSE))))
    on <- do.call(predict, args)
    off <- do.call(predict, c(args, list(fit_correct = FALSE)))
    expect_equal(as.numeric(on - off), as.numeric(weights %*% fit$laplace$delta), tolerance = 1e-9)
    link <- do.call(predict, c(args, list(se.fit = TRUE)))
    response <- do.call(predict, c(args, list(type = "response", se.fit = TRUE, delta = TRUE)))
    expect_equal(response$fit, exp(link$fit + link$se.fit^2 / 2))
    expect_equal(response$se.fit, exp(link$fit) * link$se.fit)
    set.seed(4)
    offdraws <- do.call(conditional, c(args, list(samples = 4, output = c("newdata", "beta", "latent", "object"), fit_correct = FALSE)))
    set.seed(4)
    ondraws <- do.call(conditional, c(args, list(samples = 4, output = c("newdata", "beta", "latent", "object"))))
    expect_identical(ondraws$object, offdraws$object)
    expect_equal(ondraws$latent - offdraws$latent, matrix(fit$laplace$delta, length(fit$laplace$delta), 4), ignore_attr = TRUE)
    expect_equal(ondraws$newdata - offdraws$newdata, matrix(on - off, length(on), 4), ignore_attr = TRUE)
    set.seed(4)
    expect_identical(do.call(conditional, c(oldargs, list(samples = 4, output = c("newdata", "beta", "latent", "object")))), offdraws)
    for (bad in list(NA, 1, c(TRUE, FALSE), NULL)) {
      expect_error(do.call(predict, c(args, list(fit_correct = bad))), "fit_correct")
      expect_error(do.call(conditional, c(args, list(fit_correct = bad))), "fit_correct")
    }
  }
})

test_that("response means use lognormal moments and the logit probit approximation", {
  grid <- expand.grid(m = c(-4, 0, 3), s = c(0, .3, 1, 2, 3, 10))
  for (family in c("poisson", "nbinomial", "Gamma", "inverse.gaussian", "binomial", "beta")) {
    sizes <- rep(c(1, 4, 9), length.out = nrow(grid))
    actual <- invlink_mean(grid$m, grid$s^2, family, sizes)
    if (family %in% c("binomial", "beta")) {
      reference <- mapply(function(m, s) integrate(function(z) plogis(m + s * z) * dnorm(z),
        -Inf, Inf, rel.tol = 1e-12)$value, grid$m, grid$s)
      probability <- if (family == "binomial") actual / sizes else actual
      expect_lt(max(abs(probability - reference)), .015)
      expect_identical(invlink_mean(grid$m, grid$s^2, family), plogis(grid$m / sqrt(1 + pi * grid$s^2 / 8)))
      expect_equal(actual, invlink_mean(grid$m, grid$s^2, family) * if (family == "binomial") sizes else 1)
      expect_true(all(actual >= 0 & actual <= if (family == "binomial") sizes else 1))
    } else {
      expect_identical(actual, exp(grid$m + grid$s^2 / 2))
      modest <- grid$s <= 3
      reference <- mapply(function(m, s) integrate(function(z) exp(m + s * z - z^2 / 2) / sqrt(2 * pi),
        -Inf, Inf, rel.tol = 1e-12)$value, grid$m[modest], grid$s[modest])
      expect_equal(actual[modest], reference, tolerance = 1e-7)
    }
    zero <- grid$s == 0
    expect_equal(actual[zero], invlink(grid$m[zero], family, sizes[zero]), tolerance = 1e-15)
    expect_equal(invlink_mean(grid$m, 0, family, sizes), invlink(grid$m, family, sizes), tolerance = 1e-15)
    if (family %in% c("binomial", "beta")) expect_identical(invlink_mean(grid$m, 0, family), plogis(grid$m))
  }
})

test_that("iid and local prediction retain mode curvature", {
  for (iid in c(FALSE, TRUE)) {
    fixture <- laplace_fixture(iid = iid)
    fit <- fixture$fit
    for (local in list(FALSE, list(method = "distance", size = 20))) {
      off <- predict(fit, fixture$new, local = local, se.fit = TRUE, fit_correct = FALSE)
      on <- predict(fit, fixture$new, local = local, se.fit = TRUE)
      expect_identical(on$se.fit, off$se.fit)
      expect_gt(max(abs(on$fit - off$fit)), 0)
    }
    expect_equal(augment(fit, newdata = fixture$new, fit_correct = FALSE)$.fitted,
      as.numeric(predict(fit, fixture$new, fit_correct = FALSE)))
  }
})

mean_correction_fixture <- function(random = FALSE, areal = FALSE, iid = FALSE) {
  d <- data.frame(x = seq(-1, 1, length.out = 30), cx = seq_len(30),
    y = rep(c(0, 1, 2, 1, 0, 2), 5), group = rep(letters[1:5], 6), off = seq(-.1, .2, length.out = 30))
  args <- list(formula = y ~ x + offset(off), data = d, family = "poisson")
  if (random) {
    args$random <- ~ (1 | group)
    args$randcov_initial <- randcov_initial(group = .15, known = "given")
  }
  if (areal) {
    args$data$y[28:30] <- NA
    args$W <- 1 * (abs(outer(1:30, 1:30, "-")) == 1)
    args$spcov_initial <- spcov_initial("car", de = .2, ie = .1, range = .2, known = "given")
    fit <- do.call(spgautor, args)
  } else {
    args$xcoord <- "cx"
    args$spcov_initial <- if (iid) spcov_initial("ie", ie = .2, known = "given") else
      spcov_initial("exponential", de = .2, ie = .1, range = 3, known = "given")
    fit <- do.call(spglm, args)
  }
  new <- d[1:3, ]; new$cx <- c(8.5, 15.5, 22.5)
  list(fit = fit, new = new)
}

test_that("the mean correction improves an importance-sampling conditional mean", {
  z <- mean_correction_fixture(); fit <- z$fit
  X <- model.matrix(fit); Si <- solve(as.matrix(covmatrix(fit)))
  P <- Si - Si %*% X %*% solve(crossprod(X, Si %*% X), crossprod(X, Si))
  eta <- fitted(fit, type = "link"); offset <- model.offset(model.frame(fit)); a <- eta - offset
  V <- solve(P + diag(exp(eta)))
  set.seed(439)
  n <- length(a); samples <- 100000
  noise <- matrix(rnorm(n * samples), n, samples)
  w <- a + t(chol(V)) %*% noise
  log_target <- colSums(fit$y * (w + offset) - exp(w + offset)) - .5 * colSums(w * (P %*% w))
  weights <- log_target + .5 * colSums(noise^2)
  weights <- exp(weights - max(weights)); weights <- weights / sum(weights)
  mean_reference <- as.numeric(w %*% weights)
  expect_gt(1 / sum(weights^2), 5000)
  expect_lt(mean(abs(a + fit$laplace$delta - mean_reference)), mean(abs(a - mean_reference)))
})

test_that("all GLM simulation routes shift centers without changing random errors", {
  for (random in c(FALSE, TRUE)) {
    for (areal in c(FALSE, TRUE)) {
      z <- mean_correction_fixture(random, areal)
      fit <- z$fit
      args <- if (areal) list(object = fit) else list(object = fit, newdata = z$new)
      paths <- if (areal) list(NULL) else list(FALSE,
        list(approximation = "low-rank", method_base = "all", method_new = "all"),
        list(approximation = "vecchia", size = 5, ordering = "none"))
      for (path in paths) {
        call <- c(args, list(samples = 8, output = c("newdata", "beta", "latent", "object")))
        if (!areal) call$local <- path
        if (random && (areal || identical(path, FALSE) || identical(path$approximation, "vecchia"))) {
          call$output <- c(call$output, "random")
        }
        set.seed(71); on <- do.call(conditional, call)
        set.seed(71); off <- do.call(conditional, c(call, list(fit_correct = FALSE)))
        expect_identical(on$object, off$object)
        for (component in setdiff(names(on), "object")) {
          difference <- on[[component]] - off[[component]]
          expect_equal(difference, matrix(difference[, 1], nrow(difference), ncol(difference)), ignore_attr = TRUE, tolerance = 1e-10)
        }
        expect_equal(as.numeric(on$latent[, 1] - off$latent[, 1]), unname(fit$laplace$delta), tolerance = 1e-10)
        expect_equal(as.numeric(on$beta[, 1] - off$beta[, 1]), unname(fit$laplace$betahat - coef(fit)), tolerance = 1e-10)
      }
      set.seed(871)
      draws <- do.call(conditional, c(args, list(samples = 5000)))
      prediction <- as.numeric(do.call(predict, args))
      mcse <- apply(draws, 1, sd) / sqrt(ncol(draws))
      expect_true(all(abs(rowMeans(draws) - prediction) < 5 * mcse))
    }
  }
  z <- mean_correction_fixture(iid = TRUE)
  set.seed(74); on <- conditional(z$fit, z$new, samples = 8, output = c("beta", "latent"))
  set.seed(74); off <- conditional(z$fit, z$new, samples = 8, output = c("beta", "latent"), fit_correct = FALSE)
  expect_equal(on$latent - off$latent, matrix(z$fit$laplace$delta, z$fit$n, 8), ignore_attr = TRUE)
})

cv_correction_fixture <- function(family = "poisson", areal = FALSE, partition = FALSE, random = FALSE) {
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

test_that("nonlocal GLM CV propagates the correction through held-out kriging weights", {
  for (family in c("poisson", "nbinomial", "binomial")) for (areal in c(FALSE, TRUE)) {
    fit <- cv_correction_fixture(family, areal)$fit
    X <- model.matrix(fit); S <- as.matrix(covmatrix(fit)); eta <- fitted(fit, "link")
    offset <- model.offset(model.frame(fit)); w <- eta - offset
    old <- fit; old$laplace <- NULL
    for (folds in list(seq_len(fit$n), rep(1:3, length.out = fit$n))) {
      expected <- shift <- numeric(fit$n)
      for (F in split(seq_len(fit$n), folds)) {
        Q <- solve(S[-F, -F, drop = FALSE]); Xm <- X[-F, , drop = FALSE]
        B <- solve(crossprod(Xm, Q %*% Xm), crossprod(Xm, Q))
        L <- X[F, , drop = FALSE] %*% B + S[F, -F, drop = FALSE] %*% Q %*% (diag(nrow(Xm)) - Xm %*% B)
        expected[F] <- as.numeric(L %*% (w[-F] + fit$laplace$delta[-F])) + offset[F]
        shift[F] <- as.numeric(L %*% fit$laplace$delta[-F])
      }
      args <- list(object = fit, folds_index = folds, cv_predict = TRUE, se.fit = TRUE, local = FALSE)
      on <- do.call(kcv, args)
      off <- do.call(kcv, c(args, list(fit_correct = FALSE)))
      expect_equal(on$cv_predict, expected, tolerance = 1e-9)
      expect_equal(on$cv_predict - off$cv_predict, shift, tolerance = 1e-9)
      expect_identical(on$se.fit, off$se.fit)
      args$object <- old
      expect_identical(do.call(kcv, args), off)
      response <- invlink_mean(expected, on$se.fit^2, family, fit$size)
      expect_equal(on$stats$bias, mean(fit$y - response))
      expect_equal(on$stats$MSPE, mean((fit$y - response)^2))
      expect_equal(on$stats$RMSPE, sqrt(mean((fit$y - response)^2)))
      rp <- kcv(fit, folds_index = folds, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE)
      derivative <- if (family == "binomial") fit$size * plogis(expected) * (1 - plogis(expected)) else exp(expected)
      expect_equal(rp$cv_predict, response)
      expect_equal(rp$se.fit, derivative * on$se.fit)
      if (length(unique(folds)) == fit$n) for (correction in c(FALSE, TRUE)) {
        loo <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, fit_correct = correction)
        expect_equal(loo, if (correction) on else off, tolerance = 1e-9)
        expect_identical(kcv(fit, k = fit$n, cv_predict = TRUE, se.fit = TRUE, fit_correct = correction), loo)
      }
    }
    for (bad in list(NA, 1, c(TRUE, FALSE), NULL)) {
      expect_error(loocv(fit, fit_correct = bad), "fit_correct")
      expect_error(kcv(fit, fit_correct = bad), "fit_correct")
    }
  }
})

test_that("local GLM LOO subsets corrected centers and mode curvature together", {
  fit <- cv_correction_fixture("binomial", partition = TRUE, random = TRUE)$fit
  X <- model.matrix(fit); S <- as.matrix(covmatrix(fit)); d <- fit$obdata
  w <- fitted(fit, "link") + fit$laplace$delta - d$off
  beta <- fit$laplace$betahat
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
    off <- loocv(fit, local = ctl, cv_predict = TRUE, se.fit = TRUE, fit_correct = FALSE)
    expect_equal(on$cv_predict, expected, tolerance = 1e-9)
    expect_identical(on$se.fit, off$se.fit)
    response <- loocv(fit, local = ctl, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE)
    prob <- plogis(on$cv_predict)
    expect_equal(response$cv_predict, invlink_mean(on$cv_predict, on$se.fit^2, fit$family, fit$size))
    expect_equal(response$se.fit, fit$size * prob * (1 - prob) * on$se.fit)
    old <- fit; old$laplace <- NULL
    expect_identical(loocv(old, local = ctl, cv_predict = TRUE, se.fit = TRUE), off)
  }
})

test_that("local kcv uses each training fold's correction and retains old-object fallback", {
  z <- cv_correction_fixture("binomial", partition = TRUE, random = TRUE)
  z$args$local <- list(index = rep(1:3, 8), var_adjust = "none")
  fit <- do.call(spglm, z$args)
  folds <- rep(1:3, each = 8)
  ctl <- list(method = "distance", size = 8)
  for (correction in c(FALSE, TRUE)) {
    expected <- se <- numeric(fit$n)
    for (F in split(seq_len(fit$n), folds)) {
      args <- z$args
      args$data$y[F] <- NA
      args$local$index <- fit$local_index[-F]
      train <- do.call(spglm, args)
      expect_gt(max(abs(train$laplace$delta - fit$laplace$delta[-F])), 1e-5)
      train$vcov$fixed <- fit$vcov$fixed
      pred <- predict(train, se.fit = TRUE, local = ctl, fit_correct = correction)
      expected[F] <- pred$fit; se[F] <- pred$se.fit
    }
    actual <- kcv(fit, folds_index = folds, local = ctl, cv_predict = TRUE, se.fit = TRUE, fit_correct = correction)
    expect_equal(actual$cv_predict, expected, tolerance = 1e-9)
    expect_equal(actual$se.fit, se, tolerance = 1e-9)
    if (!correction) {
      old <- fit; old$laplace <- NULL
      expect_identical(kcv(old, folds_index = folds, local = ctl, cv_predict = TRUE, se.fit = TRUE), actual)
      off <- actual
    } else expect_identical(actual$se.fit, off$se.fit)
  }
})

test_that("GLM CV corrections agree between serial and parallel paths", {
  for (areal in c(FALSE, TRUE)) {
    z <- cv_correction_fixture(areal = areal)
    fit <- z$fit
    if (!areal) {
      z$args$local <- list(index = rep(1:3, 8), var_adjust = "none")
      fit <- do.call(spglm, z$args)
    }
    methods <- if (areal) "all" else c("all", "distance")
    for (method in methods) for (correction in c(FALSE, TRUE)) {
      args <- list(object = fit, cv_predict = TRUE, se.fit = TRUE, fit_correct = correction,
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

test_that("response prediction averages the inverse link and preserves interval quantiles", {
  for (family in c("poisson", "binomial")) for (areal in c(FALSE, TRUE)) {
    z <- cv_correction_fixture(family, areal)
    fit <- z$fit
    args <- list(object = fit)
    if (!areal) {
      args$newdata <- z$args$data[1:3, ]
      args$newdata$cx <- args$newdata$cx + .02
    }
    sizes <- if (areal) z$args$data$trials[22:24] else c(2, 5, 8)
    args$newdata_size <- sizes
    methods <- if (areal) list(FALSE) else list(FALSE, list(method = "distance", size = 10))
    for (local in methods) for (correct_var in c(FALSE, TRUE)) {
      call <- c(args, list(local = local, var_correct = correct_var))
      link <- do.call(predict, c(call, list(se.fit = TRUE)))
      expected <- invlink_mean(link$fit, link$se.fit^2, family, sizes)
      for (interval in c("none", "prediction")) {
        link_bounds <- do.call(predict, c(call, list(interval = interval, se.fit = TRUE)))
        response <- do.call(predict, c(call, list(type = "response", interval = interval, se.fit = TRUE, delta = TRUE)))
        no_se <- do.call(predict, c(call, list(type = "response", interval = interval)))
        expect_equal(if (interval == "none") response$fit else response$fit[, "fit"], expected)
        expect_identical(no_se, response$fit)
        expect_equal(response$se.fit, get_delta_se(link$fit, link$se.fit, family, sizes))
        if (interval == "prediction") {
          expect_equal(response$fit[, "lwr"], invlink(link_bounds$fit[, "lwr"], family, sizes))
          expect_equal(response$fit[, "upr"], invlink(link_bounds$fit[, "upr"], family, sizes))
          expect_false("se.fit" %in% names(no_se))
        }
        old_call <- call; old_call$object$laplace <- NULL
        expect_identical(do.call(predict, c(old_call, list(type = "response", interval = interval))),
          do.call(predict, c(call, list(type = "response", interval = interval, fit_correct = FALSE))))
      }
      confidence <- do.call(predict, c(call, list(interval = "confidence", se.fit = TRUE)))
      response <- do.call(predict, c(call, list(type = "response", interval = "confidence", se.fit = TRUE)))
      expect_equal(response$fit[, "fit"], invlink(confidence$fit[, "fit"], family, sizes))
      expect_identical(response$se.fit, confidence$se.fit)
    }
    if (!areal) {
      augmented <- augment(fit, newdata = args$newdata, type.predict = "response", newdata_size = sizes)
      expect_equal(augmented$.fitted, as.numeric(do.call(predict, c(args, list(type = "response")))))
    }
  }
})

test_that("response CV computes the mean with or without requested standard errors", {
  for (family in c("poisson", "binomial")) for (areal in c(FALSE, TRUE)) {
    z <- cv_correction_fixture(family, areal)
    if (!areal) {
      z$args$local <- list(index = rep(1:3, 8), var_adjust = "none")
      z$fit <- do.call(spglm, z$args)
    }
    fit <- z$fit
    methods <- if (areal) list(FALSE) else list(FALSE, list(method = "distance", size = 10))
    for (local in methods) for (fun in list(loocv, kcv)) {
      call <- list(object = fit, local = local)
      if (identical(fun, kcv)) call$folds_index <- rep(1:3, length.out = fit$n)
      link <- do.call(fun, c(call, list(cv_predict = TRUE, se.fit = TRUE)))
      expected <- invlink_mean(link$cv_predict, link$se.fit^2, family, fit$size)
      response <- do.call(fun, c(call, list(cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE)))
      no_se <- do.call(fun, c(call, list(cv_predict = TRUE, type = "response")))
      stats <- do.call(fun, call)
      expect_equal(response$cv_predict, expected)
      expect_identical(no_se$cv_predict, response$cv_predict)
      expect_false("se.fit" %in% names(no_se))
      expect_identical(no_se$stats, response$stats)
      expect_identical(stats, response$stats)
      expect_equal(stats$bias, mean(fit$y - expected))
      expect_equal(stats$MSPE, mean((fit$y - expected)^2))
      expect_equal(stats$RMSPE, sqrt(mean((fit$y - expected)^2)))
      expect_equal(response$se.fit, get_delta_se(link$cv_predict, link$se.fit, family, fit$size))
      old_call <- call; old_call$object$laplace <- NULL
      expect_identical(do.call(fun, c(old_call, list(cv_predict = TRUE, type = "response"))),
        do.call(fun, c(call, list(cv_predict = TRUE, type = "response", fit_correct = FALSE))))
      if (family == "binomial") {
        call$local <- if (isFALSE(local)) list(method = "all") else local
        call$local$parallel <- TRUE; call$local$ncores <- 2
        expect_equal(do.call(fun, c(call, list(cv_predict = TRUE, type = "response"))), no_se, tolerance = 1e-12)
      }
    }
  }
})

test_that("conditional response draws agree within the response-mean approximation", {
  for (family in c("poisson", "binomial")) {
    z <- cv_correction_fixture(family)
    new <- z$args$data[1:3, ]; new$cx <- new$cx + .02
    sizes <- c(2, 5, 8)
    expected <- predict(z$fit, new, type = "response", newdata_size = sizes)
    set.seed(1921)
    draws <- conditional(z$fit, new, type = "response", newdata_size = sizes, samples = 20000)
    mcse <- apply(draws, 1, sd) / sqrt(ncol(draws))
    approximation <- if (family == "binomial") .015 * sizes else 0
    expect_true(all(abs(rowMeans(draws) - expected) < 4 * mcse + approximation))
  }
})
