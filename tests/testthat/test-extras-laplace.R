skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

test_that("Laplace likelihoods reproduce exact Gaussian ML and REML including constants", {
  n <- 12
  X <- cbind(1, seq(-1, 1, length.out = n))
  y <- sin(seq_len(n)) + X %*% c(.3, .7)
  offset <- seq(-.2, .2, length.out = n)^2
  logdet <- function(x) as.numeric(determinant(x, logarithm = TRUE)$modulus)
  for (path in c("point", "blocks", "areal")) {
    index <- if (path == "blocks") rep(1:3, each = 4) else rep(1, n)
    rows <- split(seq_len(n), index)
    for (scale in c(.2, 2)) {
      Sigma <- scale * (.6^abs(outer(seq_len(n), seq_len(n), "-")) + diag(.3, n))
      Sigma[outer(index, index, "!=")] <- 0
      covs <- lapply(rows, function(i) Sigma[i, i, drop = FALSE])
      precision <- Matrix::Matrix(solve(Sigma), sparse = TRUE)
      for (noise in c(.1, 1, 3)) {
        obj <- list(family = "gaussian", X = X, y = y, X_list = lapply(rows, function(i) X[i, , drop = FALSE]),
          y_list = lapply(rows, function(i) y[i, , drop = FALSE]), offset = offset,
          local_index = index, diagtol = 0, parallel = FALSE)
        V <- Sigma + diag(noise, n)
        C <- solve(crossprod(X, solve(V, X)))
        beta <- C %*% crossprod(X, solve(V, y - offset))
        residual <- y - offset - X %*% beta
        ml <- n * log(2 * pi) + logdet(V) + as.numeric(crossprod(residual, solve(V, residual)))
        for (method in c("ml", "reml")) {
          products <- with_mocked_bindings({
            fun <- if (path == "areal") laploglik_products.car else laploglik_products.exponential
            fun(NULL, noise, obj, method, NULL, NULL)
          }, get_cov_matrix_list = function(...) covs,
          spautor_cov_matrixInv = function(...) list(SigInv = precision, Sigldet = logdet(Sigma)),
          get_w_init = function(family, y, dispersion) rep(0, length(y)),
          get_d = function(family, w, y, size, dispersion) (y - w) / dispersion,
          get_D = function(family, w, y, size, dispersion) Matrix::Diagonal(n, -1 / dispersion),
          get_l00 = function(family, w, y, size, dispersion) -2 * sum(dnorm(as.numeric(y), as.numeric(w), sqrt(dispersion), log = TRUE)))
          expected <- ml + if (method == "reml") logdet(solve(C)) - ncol(X) * log(2 * pi) else 0
          expect_equal(get_minustwolaploglik(products, method, n, ncol(X)), expected,
            tolerance = 1e-10, info = paste(path, scale, noise, method))
        }
      }
    }
  }
})

test_that("independent-error ML approaches ordinary GLM likelihoods as variance shrinks", {
  set.seed(321)
  n <- 60
  dat <- data.frame(x = rnorm(n), off = runif(n, -.2, .2))
  dat$count <- rpois(n, exp(.3 + .4 * dat$x + dat$off))
  dat$success <- rbinom(n, 8, plogis(.3 + .4 * dat$x + dat$off))
  for (family in c("poisson", "binomial")) {
    formula <- if (family == "poisson") count ~ x + offset(off) else cbind(success, 8 - success) ~ x + offset(off)
    ordinary <- glm(formula, family = family, data = dat)
    gaps <- numeric()
    for (v in c(.01, .001, .0001)) {
      fit <- suppressWarnings(spglm(formula, family = as.character(family), data = dat, estmethod = "ml", local = FALSE,
        spcov_initial = spcov_initial("ie", ie = v, known = "given")))
      gaps <- c(gaps, abs(-2 * as.numeric(logLik(fit) - logLik(ordinary))))
      expect_lt(tail(gaps, 1), n * v)
      expect_equal(as.numeric(AIC(fit) - AIC(ordinary)), -2 * as.numeric(logLik(fit) - logLik(ordinary)), tolerance = 1e-9)
    }
    expect_true(all(diff(gaps) < 0))
    none <- suppressWarnings(spglm(formula, family = as.character(family), data = dat, estmethod = "ml", local = FALSE, spcov_type = "none"))
    expect_equal(as.numeric(logLik(none)), as.numeric(logLik(fit)), tolerance = 1e-10)
  }
})

test_that("Laplace determinants use the right curvature without changing fixed-parameter inference", {
  set.seed(59)
  n <- 24
  dat <- data.frame(x = rnorm(n), off = runif(n, -.2, .2), cx = runif(n), cy = runif(n),
    g = factor(rep(1:6, each = 4)), part = factor(rep(1:2, each = 12)))
  eta <- .3 + .2 * dat$x + dat$off + rnorm(n, sd = .3)
  dat$poisson <- rpois(n, exp(eta))
  dat$binomial <- rbinom(n, 8, plogis(eta))
  dat$nbinomial <- rnbinom(n, mu = exp(eta), size = 4)
  dat$Gamma <- rgamma(n, shape = 4, rate = 4 / exp(eta))
  dat$inverse.gaussian <- exp(eta) * rep(c(.7, .9, 1.1, 1.3), 6)
  dat$beta <- .25 + .5 * plogis(eta)
  W <- 1 * (abs(outer(seq_len(n), seq_len(n), "-")) == 1)
  point <- get_w_and_H_spglm
  areal <- get_w_and_H_spgautor
  original_D <- get_D
  logdet <- function(x) as.numeric(determinant(x, logarithm = TRUE)$modulus)
  for (family in c("poisson", "binomial", "nbinomial", "Gamma", "inverse.gaussian", "beta")) {
    dat$response <- dat[[family]]
    formula <- if (family == "binomial") cbind(response, 8 - response) ~ x + offset(off) else response ~ x + offset(off)
    for (path in c("point", "blocks", "random", "anisotropic", "partition", "car", "sar")) {
      args <- list(formula = formula, data = dat, family = family)
      is_areal <- path %in% c("car", "sar")
      if (is_areal) {
        args$W <- W
        args$spcov_initial <- spcov_initial(path, de = .3, ie = .2, range = .2, known = "given")
      } else {
        args$xcoord <- "cx"
        args$ycoord <- "cy"
        args$local <- if (path == "blocks") list(index = rep(1:3, each = 8)) else FALSE
        args$spcov_initial <- spcov_initial("exponential", de = .3, ie = .2, range = .3,
          rotate = if (path == "anisotropic") .3 else 0, scale = if (path == "anisotropic") .7 else 1, known = "given")
        if (path == "anisotropic") args$anisotropy <- TRUE
        if (path == "partition") args$partition_factor <- ~part
      }
      if (path == "random") {
        args$random <- ~g
        args$randcov_initial <- randcov_initial(g = .1, known = "given")
      }
      if (!family %in% c("poisson", "binomial")) args$dispersion_initial <- dispersion_initial(as.character(family), 4, known = "given")
      fits <- list()
      for (method in c("ml", "reml")) {
        args$estmethod <- method
        checked <- FALSE
        check_mode <- function(original, ...) {
          call_args <- list(...)
          D <- NULL
          call_args$ret_mHInv <- TRUE
          result <- with_mocked_bindings(do.call(original, call_args), get_D = function(...) {
            D <<- original_D(...)
            D
          })
          if (!checked) {
            checked <<- TRUE
            precision <- if (is_areal) call_args$SigInv else Matrix::bdiag(call_args$SigInv_list)
            expected <- if (method == "ml") logdet(precision - D) else -logdet(result$mHInv)
            # Preserve the existing block REML ridge convention (deferred separately).
            if (!is_areal && length(call_args$SigInv_list) > 1 && method == "reml") {
              CInv <- call_args$cov_betahat_Inv
              diag(CInv) <- diag(CInv) + call_args$data_object$diagtol
              expected <- expected + logdet(call_args$cov_betahat) + logdet(CInv)
            }
            expect_equal(result$mHldet, expected, tolerance = 1e-9, info = paste(family, path, method))
          }
          result
        }
        fits[[method]] <- with_mocked_bindings(do.call(if (is_areal) spgautor else spglm, args),
          get_w_and_H_spglm = function(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
            check_mode(point, data_object = data_object, dispersion = dispersion, SigInv_list = SigInv_list,
              SigInv_X = SigInv_X, cov_betahat = cov_betahat, cov_betahat_Inv = cov_betahat_Inv, estmethod = estmethod)
          },
          get_w_and_H_spgautor = function(data_object, dispersion, SigInv, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
            check_mode(areal, data_object = data_object, dispersion = dispersion, SigInv = SigInv,
              SigInv_X = SigInv_X, cov_betahat = cov_betahat, cov_betahat_Inv = cov_betahat_Inv, estmethod = estmethod)
          })
        expect_true(checked)
      }
      expect_identical(coef(fits$ml), coef(fits$reml))
      expect_identical(fitted(fits$ml, type = "link"), fitted(fits$reml, type = "link"))
      expect_identical(vcov(fits$ml), vcov(fits$reml))
      expect_identical(vcov(fits$ml, var_correct = FALSE), vcov(fits$reml, var_correct = FALSE))
      expect_identical(fitted(fits$ml), fitted(fits$reml))
      set.seed(92)
      ml_draws <- conditional(fits$ml, output = c("latent", "beta"), samples = 5)
      set.seed(92)
      reml_draws <- conditional(fits$reml, output = c("latent", "beta"), samples = 5)
      expect_identical(ml_draws, reml_draws)
    }
  }
})

run_ie_search <- function(args, controls, spcov) {
  optim_with_ie_bound(args, controls, ie_only_lower(args$par, spcov,
    args$data_object$diagtol, args$spcov_profiled))
}
laplace_search_fixture <- function(fn, initial = .3, de = 0, diagtol = 1e-4) {
  spcov <- spcov_initial("exponential", de = de, ie = initial, range = 1,
    known = c("de", "range"))
  args <- list(par = c(ie_log = log(initial)), fn = fn,
    data_object = list(diagtol = diagtol), spcov_profiled = FALSE)
  controls <- check_optim_method(args$par, get_optim_dotlist())
  list(args = args, controls = controls, spcov = spcov)
}

test_that("scalar ie searches use the model-derived covariance floor", {
  for (policy in list(c(0, 1e-4), c(50, 1e-4), c(2, .03))) {
    floor <- max(policy[1] / 1e4, policy[2])
    fn <- function(par, ...) (max(exp(par), floor) - .7)^2
    z <- laplace_search_fixture(fn, de = policy[1], diagtol = policy[2])
    fit <- run_ie_search(z$args, z$controls, z$spcov)
    expect_equal(unname(exp(fit$par)), .7, tolerance = 1e-5)
    expect_lt(fit$value, 1e-10)

    z$args$fn <- function(par, ...) max(exp(par), floor)
    z$args$par[] <- -Inf
    fit <- run_ie_search(z$args, z$controls, z$spcov)
    expect_equal(unname(exp(fit$par)) / floor, 1, tolerance = 1e-6)
    expect_equal(fit$value / floor, 1, tolerance = 1e-6)
    expect_identical(fit$convergence, 0L)
  }
})

test_that("bounded searches preserve optimizer controls and diagnostics", {
  z <- laplace_search_fixture(function(par, ...) -(par - 1)^2)
  z$controls$hessian <- TRUE
  z$controls$control$fnscale <- -1
  controls <- z$controls
  controls$lower <- ie_only_lower(z$args$par, z$spcov, 1e-4)
  fit <- run_ie_search(z$args, z$controls, z$spcov)
  expect_identical(fit, do.call("optim", c(z$args, controls)))

  local_mocked_bindings(optim = function(...) {
    list(par = c(ie_log = 0), value = 1, convergence = 7L,
      counts = c(functions = 12L), message = "optimizer diagnostic")
  }, .package = "stats")
  z$controls$hessian <- FALSE
  z$controls$control$fnscale <- 1
  z$args$fn <- function(par, ...) exp(par)
  fit <- run_ie_search(z$args, z$controls, z$spcov)
  expect_identical(fit$par, c(ie_log = 0))
  expect_identical(fit$value, 1)
  expect_identical(fit$convergence, 7L)
  expect_identical(fit$message, "optimizer diagnostic")
  expect_identical(fit$counts, c(functions = 12L))
})

test_that("scalar search guards leave unrelated optimization unchanged", {
  z <- laplace_search_fixture(function(par, ...) sum((par - 1)^2))
  for (name in c("range_log", "dispersion_log", "de_log")) {
    args <- z$args
    names(args$par) <- name
    expect_identical(run_ie_search(args, z$controls, z$spcov),
      do.call("optim", c(args, z$controls)))
  }
  for (type in c("car", "sar", "none")) {
    spcov <- z$spcov
    class(spcov) <- type
    expect_identical(run_ie_search(z$args, z$controls, spcov),
      do.call("optim", c(z$args, z$controls)))
  }
  z$args$par <- c(ie_log = 0, dispersion_log = 0)
  z$controls <- get_optim_dotlist()
  expect_identical(run_ie_search(z$args, z$controls, z$spcov),
    do.call("optim", c(z$args, z$controls)))
})

test_that("invalid bounds and unevaluable likelihoods still fail", {
  z <- laplace_search_fixture(function(par, ...) stop("unevaluable likelihood"))
  expect_error(run_ie_search(z$args, z$controls, z$spcov), "unevaluable likelihood")
  for (floor in c(Inf, NA_real_, exp(51))) {
    z$args$data_object$diagtol <- floor
    expect_error(run_ie_search(z$args, z$controls, z$spcov), "Invalid bounds")
  }
  z$args$data_object$diagtol <- 1e-4
  for (upper in list(NULL, numeric(), c(0, 1), Inf, NA_real_, log(1e-4))) {
    z$controls$upper <- upper
    expect_error(run_ie_search(z$args, z$controls, z$spcov), "Invalid bounds")
  }
  expect_error(optim_with_ie_bound(list(par = c(range_log = 0)), list(upper = 50), -10), "Invalid bounds")
})

test_that("profiled ie searches use the proportion floor", {
  z <- laplace_search_fixture(function(par, ...) (plogis(par) - .4)^2, diagtol = 0)
  z$args$par <- c(ie_prop_logodds = qlogis(.3))
  z$args$spcov_profiled <- TRUE
  lower <- ie_only_lower(z$args$par, z$spcov, 0, TRUE)
  expect_equal(lower, log(1e-4))
  fit <- run_ie_search(z$args, z$controls, z$spcov)
  expect_equal(unname(plogis(fit$par)), .4, tolerance = 1e-6)
  z$args$fn <- function(par, ...) max(plogis(par), plogis(lower))
  fit <- run_ie_search(z$args, z$controls, z$spcov)
  expect_equal(unname(exp(fit$par - lower)), 1, tolerance = 1e-6)
  expect_null(ie_only_lower(z$args$par, z$spcov, .01, TRUE))

  z <- laplace_search_fixture(function(par, ...) (par - 1)^2, diagtol = 0)
  expect_null(ie_only_lower(z$args$par, z$spcov, 0))
  expect_identical(run_ie_search(z$args, z$controls, z$spcov),
    do.call("optim", c(z$args, z$controls)))
  expect_null(ie_only_lower(z$args$par, z$spcov, 0, TRUE))
})

test_that("estimated ie recovers the independently bracketed fixed-variance optimum", {
  set.seed(12)
  n <- 300
  dat <- data.frame(x = rnorm(n))
  dat$y <- rpois(n, exp(.5 + .3 * dat$x + rnorm(n, sd = .8)))
  for (method in c("reml", "ml")) {
    fixed <- function(theta) {
      spglm(y ~ x, family = "poisson", data = dat, estmethod = method,
        spcov_initial = spcov_initial("ie", ie = exp(theta), known = "given"), local = FALSE)
    }
    reference <- optimize(function(theta) -2 * as.numeric(logLik(fixed(theta))),
      interval = log(c(.01, 10)), tol = 1e-6)
    boundary <- suppressWarnings(fixed(log(1e-4)))
    for (start in c(1e-8, .3, 10)) {
      fit <- spglm(y ~ x, family = "poisson", data = dat, estmethod = method,
        spcov_initial = spcov_initial("ie", ie = start), local = FALSE)
      expect_equal(unname(coef(fit, type = "spcov")[["ie"]]), exp(reference$minimum), tolerance = 1e-4)
      expect_equal(fit$optim$value, reference$objective, tolerance = 1e-7)
      expect_equal(fit$optim$value, -2 * as.numeric(logLik(fit)), tolerance = 1e-10)
      expect_lt(fit$optim$value, -2 * as.numeric(logLik(boundary)) - 100)
      expect_equal(coef(fit), coef(fixed(reference$minimum)), tolerance = 1e-5)
      expect_false(fit$is_known$spcov[["ie"]])
    }
  }
})

test_that("scalar ie fits agree with fixed-parameter likelihoods across supported contexts", {
  set.seed(27)
  n <- 60
  dat <- data.frame(x = rnorm(n), off = runif(n, -.2, .2), cx = runif(n), cy = runif(n),
    g = rep(letters[1:5], length.out = n), part = factor(rep(1:2, length.out = n)))
  eta <- .5 + .25 * dat$x + dat$off + rnorm(n, sd = .8)
  dat$y <- rpois(n, exp(eta))
  dat$gam <- rgamma(n, shape = 4, rate = 4 / exp(eta))
  dat$bin <- rbinom(n, 10, plogis(eta))
  dat$nb <- rnbinom(n, size = 4, mu = exp(eta))
  dat$ig <- statmod::rinvgauss(n, mean = exp(eta), shape = 4)
  mu <- plogis(eta)
  dat$beta_extreme <- rbeta(n, mu * 4, (1 - mu) * 4)
  dat$beta <- .2 + .6 * dat$beta_extreme
  base <- list(formula = y ~ x + offset(off), family = "poisson", data = dat,
    spcov_initial = spcov_initial("ie", ie = .3), local = FALSE)
  cases <- list(poisson = base, blocks = modifyList(base, list(local = list(index = rep(1:3, each = 20)))),
    partitions = c(base, list(partition_factor = ~part)),
    random = c(base, list(random = ~g, randcov_initial = randcov_initial(g = .1, known = "given"))))
  for (family in c("Gamma", "binomial", "nbinomial", "inverse.gaussian", "beta")) {
    args <- base
    args$family <- family
    args$formula <- switch(family, Gamma = gam ~ x + offset(off),
      binomial = cbind(bin, 10 - bin) ~ x + offset(off), nbinomial = nb ~ x + offset(off),
      inverse.gaussian = ig ~ x + offset(off), beta = beta ~ x + offset(off))
    if (family != "binomial") args$dispersion_initial <- dispersion_initial(as.character(family), dispersion = 4, known = "given")
    cases[[family]] <- args
  }
  for (anis in c(FALSE, TRUE)) {
    args <- base
    args$spcov_initial <- spcov_initial("exponential", de = .1, ie = .3, range = .3,
      rotate = if (anis) .3 else 0, scale = if (anis) .7 else 1,
      known = c("de", "range", "rotate", "scale"))
    args$xcoord <- "cx"
    args$ycoord <- "cy"
    args$anisotropy <- anis
    cases[[if (anis) "anisotropic" else "spatial"]] <- args
  }
  for (label in names(cases)) {
    for (method in c("ml", "reml")) {
      args <- cases[[label]]
      args$estmethod <- method
      fit <- suppressWarnings(do.call(spglm, args))
      ie <- coef(fit, type = "spcov")[["ie"]]
      fixed <- function(variance) {
        fixed_args <- args
        fixed_args$spcov_initial$initial[["ie"]] <- variance
        fixed_args$spcov_initial$is_known[["ie"]] <- TRUE
        ref <- suppressWarnings(do.call(spglm, fixed_args))
        if (label == "anisotropic") {
          fixed_args$spcov_initial$initial[["rotate"]] <- pi - fixed_args$spcov_initial$initial[["rotate"]]
          alternate <- suppressWarnings(do.call(spglm, fixed_args))
          if (as.numeric(logLik(alternate)) > as.numeric(logLik(ref))) ref <- alternate
        }
        ref
      }
      ref <- fixed(ie)
      info <- paste(label, method)
      expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-8, info = info)
      expect_equal(coef(fit), coef(ref), tolerance = 1e-8, info = info)
      expect_equal(fit$optim$value, -2 * as.numeric(logLik(fit)), tolerance = 1e-8, info = info)
      expect_false(fit$is_known$spcov[["ie"]], info = info)
      if (label %in% c("spatial", "anisotropic")) {
        expect_identical(fit$anisotropy, label == "anisotropic")
      }
      for (candidate in c(1e-4, .3, ie / 1.05, ie * 1.05)) {
        expect_true(fit$optim$value <= -2 * as.numeric(logLik(fixed(candidate))) + 1e-5, info = info)
      }
    }
  }
  args <- cases$beta
  args$formula <- beta_extreme ~ x + offset(off)
  expect_error(suppressWarnings(do.call(spglm, args)), "Convergence problem")
  expect_error(with_mocked_bindings(suppressWarnings(do.call(spglm, args)),
    optim_with_ie_bound = function(optim_args, optim_dotlist, lower) {
      do.call("optim", c(optim_args, optim_dotlist))
    }), "Convergence problem")
})

test_that("boundary fits approach the floor and agree at the returned variance", {
  dat <- data.frame(y = rep(c(1, 2), 30))
  for (method in c("reml", "ml")) {
    fit <- suppressWarnings(spglm(y ~ 1, data = dat, family = "poisson",
      spcov_type = "ie", estmethod = method, local = FALSE))
    floor <- suppressWarnings(spglm(y ~ 1, data = dat, family = "poisson",
      spcov_initial = spcov_initial("ie", ie = 1e-4, known = "given"), estmethod = method, local = FALSE))
    ie <- coef(fit, type = "spcov")[["ie"]]
    expect_equal(ie / 1e-4, 1, tolerance = 1e-6)
    expect_equal(fit$optim$value, -2 * as.numeric(logLik(floor)), tolerance = 1e-7)
    ref <- suppressWarnings(spglm(y ~ 1, data = dat, family = "poisson",
      spcov_initial = spcov_initial("ie", ie = ie, known = "given"), estmethod = method, local = FALSE))
    expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-12)
    expect_equal(coef(fit), coef(ref), tolerance = 1e-10)
  }
})

test_that("model-specific floors agree with fixed-parameter reference searches", {
  set.seed(51)
  dat <- data.frame(y = rep(c(1, 2), 20), x = rnorm(40), cx = runif(40), cy = runif(40))
  original <- get_data_object_spglm
  for (policy in list(c(50, 1e-4), c(.1, .01), c(.1, 1e-6))) {
    floor <- max(policy[1] / 1e4, policy[2])
    fit_at <- function(ie, known = c("de", "range")) {
      with_mocked_bindings(
        suppressWarnings(spglm(y ~ x, data = dat, family = "poisson", xcoord = "cx", ycoord = "cy",
          spcov_initial = spcov_initial("exponential", de = policy[1], ie = ie, range = .1, known = known), local = FALSE)),
        get_data_object_spglm = function(...) {
          obj <- original(...)
          obj$diagtol <- policy[2]
          obj
        })
    }
    fixed <- function(theta) -2 * as.numeric(logLik(fit_at(exp(theta), "given")))
    reference <- optimize(fixed, interval = c(log(floor), log(10)), tol = 1e-6)
    target <- min(reference$objective, fixed(log(floor)))
    fit <- fit_at(.3)
    expect_equal(fit$optim$value, target, tolerance = 1e-7)
    expect_gte(coef(fit, type = "spcov")[["ie"]], floor)
    expect_equal(fit$optim$value, fixed(log(coef(fit, type = "spcov")[["ie"]])), tolerance = 1e-9)
    expect_true(fit$is_known$spcov[["de"]])
    expect_false(fit$is_known$spcov[["ie"]])
  }
})

test_that("fixed, other-scalar and multi-parameter fits retain their optimizer behavior", {
  set.seed(72)
  n <- 40
  dat <- data.frame(x = rnorm(n), cx = runif(n), cy = runif(n), g = rep(letters[1:4], each = 10))
  dat$y <- rpois(n, exp(.5 + .3 * dat$x + rnorm(n, sd = .5)))
  dat$gam <- rgamma(n, shape = 4, rate = 4 / exp(.5 + .3 * dat$x))
  base <- list(formula = y ~ x, data = dat, family = "poisson", local = FALSE,
    spcov_initial = spcov_initial("ie", ie = .3, known = "given"))
  cases <- list(fixed = base)
  cases$none <- base
  cases$none$spcov_initial <- spcov_initial("none")
  cases$dispersion <- base
  cases$dispersion$formula <- gam ~ x
  cases$dispersion$family <- "Gamma"
  cases$dispersion$dispersion_initial <- dispersion_initial("Gamma", dispersion = 4)
  cases$multiple <- cases$dispersion
  cases$multiple$spcov_initial <- spcov_initial("ie", ie = .3)
  cases$random <- c(base, list(random = ~g, randcov_initial = randcov_initial(g = .1)))
  cases$random$spcov_initial <- spcov_initial("ie", ie = .3)
  cases$range <- c(base, list(xcoord = "cx", ycoord = "cy"))
  cases$range$spcov_initial <- spcov_initial("exponential", de = .1, ie = .3, range = .3,
    known = c("de", "ie"))
  for (label in names(cases)) {
    args <- cases[[label]]
    fit <- suppressWarnings(do.call(spglm, args))
    reference <- with_mocked_bindings(suppressWarnings(do.call(spglm, args)),
      optim_with_ie_bound = function(optim_args, optim_dotlist, lower) {
        do.call("optim", c(optim_args, optim_dotlist))
      })
    for (field in c("coefficients", "optim", "is_known", "fitted", "vcov")) {
      expect_identical(fit[[field]], reference[[field]], info = paste(label, field))
    }
  }
  gaussian <- with_mocked_bindings(splm(y ~ x, data = dat, spcov_type = "ie", local = FALSE),
    optim_with_ie_bound = function(...) stop("Scalar optimizer reached"))
  expect_s3_class(gaussian, "splm")
})
