skip_on_cran()
skip_if_not(identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"))

gaussian_search_data <- function(n = 200) {
  set.seed(4)
  dat <- data.frame(cx = runif(n), cy = runif(n), x = rnorm(n))
  S <- exp(-as.matrix(dist(dat[, c("cx", "cy")])) / .2)
  dat$y <- as.vector(1 + .5 * dat$x + t(chol(S)) %*% rnorm(n) + rnorm(n, sd = sqrt(.5)))
  dat
}

test_that("Gaussian scalar searches recover fixed-parameter reference optima", {
  dat <- gaussian_search_data()
  for (method in c("reml", "ml")) {
    fixed <- function(de, ie) splm(y ~ x, data = dat, xcoord = cx, ycoord = cy,
      estmethod = method, local = FALSE,
      spcov_initial = spcov_initial("exponential", de = de, ie = ie, range = .2, known = "given"))
    for (profiled in c(FALSE, TRUE)) {
      objective <- function(theta) {
        if (!profiled) return(-2 * as.numeric(logLik(fixed(1, exp(theta)))))
        q <- plogis(theta)
        optimize(function(log_total) -2 * as.numeric(logLik(
          fixed(exp(log_total) * (1 - q), exp(log_total) * q))),
          interval = log(c(.01, 100)), tol = 1e-8)$objective
      }
      reference <- optimize(objective, interval = if (profiled) qlogis(c(.01, .99)) else log(c(.01, 10)), tol = 1e-7)
      for (start in c(1e-6, .3, 5)) {
        fit <- splm(y ~ x, data = dat, xcoord = cx, ycoord = cy, estmethod = method,
          local = FALSE, spcov_initial = spcov_initial("exponential", de = 1, ie = start,
            range = .2, known = if (profiled) "range" else c("de", "range")))
        sp <- coef(fit, type = "spcov")
        theta <- if (profiled) log(sp[["ie"]] / sp[["de"]]) else log(sp[["ie"]])
        ref <- fixed(sp[["de"]], sp[["ie"]])
        expect_equal(theta, reference$minimum, tolerance = 1e-4)
        expect_equal(fit$optim$value, reference$objective, tolerance = 1e-8)
        expect_lt(fit$optim$value, objective(log(1e-4)) - 30)
        expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-10)
        expect_equal(fit$optim$value, -2 * as.numeric(logLik(fit)), tolerance = 1e-10)
        expect_equal(coef(fit), coef(ref), tolerance = 1e-10)
        expect_equal(vcov(fit), vcov(ref), tolerance = 1e-10)
        expect_false(fit$is_known$spcov[["ie"]])
        expect_identical(fit$is_known$spcov[["de"]], !profiled)
      }
    }
  }
})

test_that("Gaussian scalar searches preserve reconstruction across fitting contexts", {
  dat <- gaussian_search_data(48)
  dat$part <- factor(rep(1:2, each = 24))
  dat$g <- rep(letters[1:6], each = 8)
  dat$off <- seq(-.3, .3, length.out = nrow(dat))
  base <- list(formula = y ~ x + offset(off), data = dat, xcoord = "cx", ycoord = "cy", local = FALSE)
  cases <- list(spatial = base,
    anisotropic = c(base, list(anisotropy = TRUE)),
    partition = c(base, list(partition_factor = ~part)),
    local = modifyList(base, list(local = list(index = rep(1:3, each = 16)))),
    random = c(base, list(random = ~g, randcov_initial = randcov_initial(g = .1, known = "given"))))
  for (label in names(cases)) for (method in c("ml", "reml")) for (profiled in c(FALSE, TRUE)) {
    if (profiled && label == "random") next
    args <- cases[[label]]
    args$estmethod <- method
    args$spcov_initial <- spcov_initial("exponential", de = 1, ie = .3, range = .2,
      rotate = if (label == "anisotropic") .3 else 0, scale = if (label == "anisotropic") .7 else 1,
      known = c(if (!profiled) "de", "range", "rotate", "scale"))
    fit <- do.call(splm, args)
    sp <- coef(fit, type = "spcov")
    fixed <- function(theta) {
      fixed_args <- args
      fixed_args$spcov_initial$initial <- sp
      fixed_args$spcov_initial$is_known[] <- TRUE
      q <- plogis(theta)
      at_total <- function(total) {
        fixed_args$spcov_initial$initial[["de"]] <- if (profiled) total * (1 - q) else sp[["de"]]
        fixed_args$spcov_initial$initial[["ie"]] <- if (profiled) total * q else exp(theta)
        ref <- do.call(splm, fixed_args)
        if (label == "anisotropic") {
          fixed_args$spcov_initial$initial[["rotate"]] <- pi - sp[["rotate"]]
          alternate <- do.call(splm, fixed_args)
          if (as.numeric(logLik(alternate)) > as.numeric(logLik(ref))) ref <- alternate
        }
        ref
      }
      if (profiled) {
        total <- optimize(function(log_total) -2 * as.numeric(logLik(at_total(exp(log_total)))),
          log(c(.001, 100)), tol = 1e-8)$minimum
        at_total(exp(total))
      } else at_total(1)
    }
    theta <- if (profiled) log(sp[["ie"]] / sp[["de"]]) else log(sp[["ie"]])
    ref <- fixed(theta)
    info <- paste(label, method, profiled)
    expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-8, info = info)
    expect_equal(coef(fit), coef(ref), tolerance = 1e-7, info = info)
    expect_equal(vcov(fit), vcov(ref), tolerance = 1e-7, info = info)
    expect_equal(fit$optim$value, -2 * as.numeric(logLik(fit)), tolerance = 1e-9, info = info)
    for (candidate in theta + c(-.05, .05)) {
      expect_lte(fit$optim$value, -2 * as.numeric(logLik(fixed(candidate))) + 1e-5)
    }
  }
})

test_that("Gaussian floor optima agree with fixed-parameter fits", {
  dat <- gaussian_search_data(40)
  dat$y <- sin(dat$cx) + cos(dat$cy)
  for (method in c("ml", "reml")) for (profiled in c(FALSE, TRUE)) {
    fit <- splm(y ~ 1, data = dat, xcoord = cx, ycoord = cy, estmethod = method,
      spcov_initial = spcov_initial("exponential", de = 1, ie = .3, range = .2,
        known = if (profiled) "range" else c("de", "range")), local = FALSE)
    sp <- coef(fit, type = "spcov")
    expect_equal(sp[["ie"]] / sp[["de"]] / 1e-4, 1, tolerance = 1e-6)
    ref <- splm(y ~ 1, data = dat, xcoord = cx, ycoord = cy, estmethod = method,
      spcov_initial = spcov_initial("exponential", de = sp[["de"]], ie = sp[["ie"]], range = .2,
        known = "given"), local = FALSE)
    expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-10)
    expect_equal(coef(fit), coef(ref), tolerance = 1e-10)
  }
})

test_that("unrelated Gaussian fits retain their optimizer behavior", {
  dat <- gaussian_search_data(40)
  dat$g <- rep(letters[1:4], each = 10)
  cases <- list(
    fixed = spcov_initial("exponential", de = 1, ie = .3, range = .2, known = "given"),
    iid = spcov_initial("ie"),
    multiple = spcov_initial("exponential", de = 1, ie = .3, range = .2),
    range = spcov_initial("exponential", de = 1, ie = .3, range = .2, known = c("de", "ie")),
    zero_floor = spcov_initial("exponential", de = 0, ie = .3, range = .2, known = c("de", "range")))
  for (label in names(cases)) {
    args <- list(formula = y ~ x, data = dat, xcoord = "cx", ycoord = "cy", local = FALSE,
      spcov_initial = cases[[label]])
    if (label == "zero_floor") {
      args$random <- ~g
      args$randcov_initial <- randcov_initial(g = .1, known = "given")
    }
    fit <- do.call(splm, args)
    ref <- with_mocked_bindings(do.call(splm, args),
      optim_with_ie_bound = function(optim_args, optim_dotlist, lower) do.call("optim", c(optim_args, optim_dotlist)))
    for (field in c("coefficients", "optim", "is_known", "fitted", "vcov")) {
      expect_identical(fit[[field]], ref[[field]], info = paste(label, field))
    }
  }
  for (type in c("car", "sar")) {
    W <- 1 * (as.matrix(dist(dat[, c("cx", "cy")])) < .4)
    diag(W) <- 0
    args <- list(formula = y ~ x, data = dat, W = W,
      spcov_initial = spcov_initial(type, de = 1, ie = .3, range = .2, known = c("de", "range")))
    fit <- do.call(spautor, args)
    ref <- with_mocked_bindings(do.call(spautor, args),
      optim_with_ie_bound = function(optim_args, optim_dotlist, lower) do.call("optim", c(optim_args, optim_dotlist)))
    expect_identical(fit$optim, ref$optim)
    expect_identical(coef(fit, type = "spcov"), coef(ref, type = "spcov"))
  }
})
