skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

test_that("prediction-grid methods do not partially match the observed method", {
  set.seed(2)
  d <- data.frame(x = runif(24), y = runif(24), z = rnorm(24))
  fit <- splm(z ~ x, d, xcoord = x, ycoord = y, ddf = "asymptotic",
    spcov_initial = spcov_initial("exponential", de = 0.8, ie = 0.2, range = 0.15, known = "given"))
  nd <- data.frame(x = runif(9), y = runif(9))
  quantities <- get_block_quantities
  calls <- list()
  testthat::local_mocked_bindings(get_block_quantities = function(object, newdata, nodes, ...) {
    calls[[length(calls) + 1L]] <<- list(rows = rownames(newdata), nodes = nodes)
    quantities(object, newdata, nodes, ...)
  }, .package = "spmodel")
  old <- options(warnPartialMatchDollar = TRUE)
  on.exit(options(old), add = TRUE)
  for (method in c("subset", "basis")) {
    settings <- list(method_new = method, size_new = 3, size = 6, ordering = "none")
    calls <- list()
    set.seed(2)
    expect_warning(a <- predict(fit, nd, block = TRUE, se.fit = TRUE, local = settings), NA)
    selected <- calls
    calls <- list()
    settings$method <- "covariance"
    set.seed(2)
    expect_warning(b <- predict(fit, nd, block = TRUE, se.fit = TRUE, local = settings), NA)
    expect_equal(a, b, tolerance = 0)
    expect_identical(selected, calls)
    expect_identical(selected[[1L]]$nodes, 1:3)
    expect_length(selected[[1L]]$rows, if (method == "subset") 3L else 9L)
  }
})

# Scaling predict(block = TRUE) to large prediction grids. The block point
# prediction stays exact on every path; only the block standard error is
# approximated when method_new = "basis" (default) with size_new < nrow(newdata),
# or when method_new = "subset".

# fixed-parameter splm so each fit is a single GLS solve
bp_model <- function(n_obs, de, ie, range, formula = z ~ 1, seed = 1) {
  set.seed(seed)
  d <- data.frame(xco = runif(n_obs), yco = runif(n_obs))
  d$s <- as.numeric(sprnorm(
    spcov_params("exponential", de = de, ie = ie, range = range),
    data = d, xcoord = xco, ycoord = yco
  ))
  # a covariate with real fine-scale spatial structure (checks x0 is not subset)
  d$x <- sin(6 * d$xco) * cos(5 * d$yco)
  d$z <- 1 + 0.5 * d$x + d$s
  splm(formula, d,
    spcov_initial = spcov_initial("exponential", de = de, ie = ie, range = range, known = "given"),
    xcoord = xco, ycoord = yco
  )
}
bp_grid <- function(G, lims = c(0, 1), seed = 99) {
  set.seed(seed)
  data.frame(
    xco = runif(G, lims[1], lims[2]), yco = runif(G, lims[1], lims[2]),
    x = NA_real_
  ) |> within(x <- sin(6 * xco) * cos(5 * yco))
}

test_that("get_block_quantities() matches the pre-change dense c0 / row-by-row s0", {
  m <- bp_model(300, de = 0.8, ie = 0.2, range = 0.15)
  gd <- bp_grid(2500)
  ie <- m$coefficients$spcov[["ie"]]

  # reference: exact c0 (dense colMeans) and exact s0 (dense pred.pred mean)
  c0_ref <- colMeans(covmatrix(m, newdata = gd, cov_type = "pred.obs"))
  s0_ref <- mean(covmatrix(m, newdata = gd, cov_type = "pred.pred"))

  bq <- spmodel:::get_block_quantities(m, gd, nodes = seq_len(nrow(gd)))
  expect_equal(bq$c0, unname(c0_ref), tolerance = 1e-10)
  expect_equal(bq$s0, s0_ref, tolerance = 1e-10)

  # row-by-row reconstruction of the removed get_bk_var()
  s0_rowwise <- mean(vapply(seq_len(nrow(gd)), function(i) {
    v <- as.vector(covmatrix(`$<-`(m, "obdata", gd), newdata = gd[i, , drop = FALSE], cov_type = "obs.pred"))
    v[i] <- v[i] + ie
    mean(v)
  }, numeric(1)))
  expect_equal(bq$s0, s0_rowwise, tolerance = 1e-10)
})

test_that("block covariance includes the stabilized ie on every grid path", {
  dat <- data.frame(
    xco = rep(seq_len(6), 4), yco = rep(seq_len(4), each = 6),
    x = sin(seq_len(24)), grp = factor(rep(c("a", "b"), 12))
  )
  dat$z <- 1 + dat$x + cos(seq_len(24))
  gd <- data.frame(
    xco = seq(6.2, 8.3, length.out = 8), yco = rep(c(1.3, 3.7), 4),
    x = cos(seq_len(8)), grp = factor(rep(c("a", "b"), 4))
  )
  G <- nrow(gd)
  nodes <- get_decorrelate_order("maxmin", gd$xco, gd$yco)$order[seq_len(4)]

  for (structure in c("plain", "partition", "random")) {
    for (ie in c(0, 1e-9, 0.2)) {
      m <- splm(z ~ x, dat, xcoord = xco, ycoord = yco, ddf = "asymptotic",
        spcov_initial = spcov_initial("exponential",
          de = 0.4, ie = ie, range = 2, known = "given"),
        partition_factor = if (structure == "partition") ~grp else NULL,
        random = if (structure == "random") ~grp else NULL,
        randcov_initial = if (structure == "random") randcov_initial(grp = 0.1, known = "given") else NULL
      )
      for (diagtol in c(0, 0.003)) {
        # Exercise an absolute floor as well as the relative floor.
        m$diagtol <- diagtol
        effective_ie <- max(ie, 4e-5, diagtol)
        V <- covmatrix(m)
        C <- covmatrix(m, gd, cov_type = "pred.obs")
        T <- covmatrix(m, gd, cov_type = "pred.pred")
        expect_equal(unname(diag(T)), rep(0.4 + effective_ie +
          if (structure == "random") 0.1 else 0, G))
        X <- model.matrix(m)
        x0 <- colMeans(model.matrix(delete.response(terms(m)), gd))
        B <- vcov(m)

        for (method in c("exact", "basis", "subset")) {
          if (method == "exact") {
            c0 <- colMeans(C)
            s0 <- mean(T)
            local <- FALSE
            bq <- get_block_quantities(m, gd, seq_len(G), chunk = 3L)
          } else if (method == "basis") {
            c0 <- colMeans(C)
            s0 <- mean(T[, nodes, drop = FALSE])
            local <- list(method = "all", method_new = method, size_new = length(nodes), ordering = "maxmin")
            bq <- get_block_quantities(m, gd, nodes, chunk = 3L)
          } else {
            c0 <- colMeans(C[nodes, , drop = FALSE])
            s0 <- mean(T[nodes, nodes])
            local <- list(method = "all", method_new = method, size_new = length(nodes), ordering = "maxmin")
            bq <- get_block_quantities(m, gd[nodes, , drop = FALSE], seq_along(nodes), chunk = 3L)
          }
          expect_equal(bq$s0, s0, tolerance = 1e-12)
          expect_equal(bq$c0, unname(c0), tolerance = 1e-12)
          if (method == "subset") {
            s0 <- s0 + effective_ie * (1 / G - 1 / length(nodes))
          }
          H <- x0 - c0 %*% solve(V, X)
          variance <- as.numeric(s0 - c0 %*% solve(V, c0) + H %*% B %*% t(H))
          got <- predict(m, gd, block = TRUE, se.fit = TRUE, local = local)
          expect_equal(as.numeric(got$se.fit^2), variance, tolerance = 1e-11)
        }
        expect_equal(coef(m, "spcov")[["ie"]], ie)
      }
    }
  }
})

test_that("chunked accumulation is independent of the chunk size", {
  m <- bp_model(200, de = 0.7, ie = 0.3, range = 0.2)
  gd <- bp_grid(1800)
  nodes <- seq_len(nrow(gd))
  base <- spmodel:::get_block_quantities(m, gd, nodes = nodes, chunk = 1000L)
  for (ch in c(1L, 37L, 250L, 5000L)) {
    alt <- spmodel:::get_block_quantities(m, gd, nodes = nodes, chunk = ch)
    expect_equal(alt$c0, base$c0, tolerance = 1e-11)
    expect_equal(alt$s0, base$s0, tolerance = 1e-11)
  }
})

test_that("method_new = 'basis' keeps the block point prediction exact", {
  m <- bp_model(400, de = 0.8, ie = 0.2, range = 0.15, formula = z ~ x)
  gd <- bp_grid(5000)
  ex <- predict(m, gd, block = TRUE, se.fit = TRUE, local = FALSE)
  for (ord in c("maxmin", "grts")) {
    set.seed(5)
    ap <- predict(m, gd,
      block = TRUE, se.fit = TRUE,
      local = list(method = "all", method_new = "basis", size_new = 1500, ordering = ord)
    )
    expect_equal(ap$fit, ex$fit, tolerance = 1e-8) # only s0 is approximate
  }
  # prediction interval endpoints also driven by the exact fit
  api <- predict(m, gd, block = TRUE, interval = "prediction", local = list(method = "all", method_new = "basis", size_new = 1500))
  exi <- predict(m, gd, block = TRUE, interval = "prediction", local = FALSE)
  expect_equal(unname(api[1, "fit"]), unname(exi[1, "fit"]), tolerance = 1e-8)
})

test_that("size_new >= nrow(newdata) reproduces the exact standard error", {
  m <- bp_model(250, de = 0.8, ie = 0.2, range = 0.2, formula = z ~ x)
  gd <- bp_grid(600)
  ex <- predict(m, gd, block = TRUE, se.fit = TRUE, local = FALSE)
  b_big <- predict(m, gd, block = TRUE, se.fit = TRUE, local = list(method = "all", method_new = "basis", size_new = 4000))
  expect_identical(b_big$fit, ex$fit)
  expect_identical(b_big$se.fit, ex$se.fit)
})

test_that("check s0 (method_new = 'basis') is accurate", {
  gd <- bp_grid(5000)
  for (rng in c(0.15, 0.40)) {
    m <- bp_model(180, de = 0.8, ie = 0.2, range = rng, formula = z ~ x)
    ex <- predict(m, gd, block = TRUE, se.fit = TRUE, local = FALSE)
    set.seed(2)
    mm <- predict(m, gd, block = TRUE, se.fit = TRUE, local = list(method = "all", method_new = "basis", size_new = 2000, ordering = "maxmin"))
    set.seed(2)
    expect_equal(mm$se.fit, ex$se.fit, tolerance = 0.1) # less accurate than grts in this example
    gr <- predict(m, gd, block = TRUE, se.fit = TRUE, local = list(method = "all", method_new = "basis", size_new = 2000, ordering = "grts"))
    expect_equal(gr$se.fit, ex$se.fit, tolerance = 0.05)
  }
})

test_that("method_new = 'subset' uses the node covariance and full-grid nugget", {
  gd <- bp_grid(6000, lims = c(0.1, 0.4)) # sub-region block
  m <- bp_model(220, de = 0.8, ie = 0.2, range = 0.15, formula = z ~ x)
  ex <- predict(m, gd, block = TRUE, se.fit = TRUE, local = FALSE)
  set.seed(4)
  sm <- predict(m, gd, block = TRUE, se.fit = TRUE, local = list(method = "all", method_new = "subset", size_new = 2000, ordering = "grts"))
  expect_equal(sm$fit, ex$fit, tolerance = 5e-3) # small block-mean error
  set.seed(4)
  nodes <- gd[get_decorrelate_order("grts", gd$xco, gd$yco)$order[seq_len(2000)], ]
  c0 <- colMeans(covmatrix(m, nodes, cov_type = "pred.obs"))
  ie <- spcov_ie_stabilized(coef(m, "spcov"), m$diagtol)
  s0 <- mean(covmatrix(m, nodes, cov_type = "pred.pred")) + ie * (1 / nrow(gd) - 1 / nrow(nodes))
  w <- solve(covmatrix(m), c0)
  h <- colMeans(model.matrix(~x, gd)) - drop(crossprod(model.matrix(m), w))
  expect_equal(sm$se.fit^2, s0 - sum(w * c0) + drop(crossprod(h, vcov(m) %*% h)),
    tolerance = 1e-10)
})

test_that("the unset-local auto-flag fires per side, only for interval none/prediction", {
  m <- bp_model(250, de = 0.8, ie = 0.2, range = 0.15)

  # small grid, small model: no approximation, no message
  gd_small <- bp_grid(400)
  expect_no_message(p_small <- predict(m, gd_small, block = TRUE, se.fit = TRUE))
  p_exact <- predict(m, gd_small, block = TRUE, se.fit = TRUE, local = FALSE)
  expect_identical(p_small$se.fit, p_exact$se.fit)

  # large grid: message names the prediction side, fit stays exact
  gd_big <- bp_grid(10200)
  expect_message(
    p_big <- predict(m, gd_big, block = TRUE, se.fit = TRUE),
    "number of prediction locations in newdata exceeds 10,000"
  )
  fit_exact_big <- predict(m, gd_big, block = TRUE, se.fit = TRUE, local = FALSE)$fit
  expect_equal(predict(m, gd_big, block = TRUE, local = FALSE)[1], p_big$fit, tolerance = 1e-8)

  # confidence intervals use neither c0 nor s0 -> never flaged
  expect_no_message(predict(m, gd_big, block = TRUE, interval = "confidence"))

  # a large observed sample flags the observed-side message (s0 still exact
  # because size_new stays Inf); patch object$n to avoid an overly large (>10k) fit
  m_bign <- m
  m_bign$n <- 10001L
  expect_message(
    predict(m_bign, gd_small, block = TRUE, se.fit = TRUE),
    "fitted model sample size exceeds 10,000"
  )

  # explicit local = FALSE always suppresses the flag
  expect_no_message(predict(m, gd_big, block = TRUE, se.fit = TRUE, local = FALSE))
})

test_that("block ordering checks: determinism, seeding, validation, coincident coords", {
  m <- bp_model(200, de = 0.8, ie = 0.2, range = 0.2)
  gd <- bp_grid(3000)
  cfg <- function(ord) list(method = "all", method_new = "basis", size_new = 800, ordering = ord)

  # maxmin is deterministic across calls
  a <- predict(m, gd, block = TRUE, se.fit = TRUE, local = cfg("maxmin"))$se.fit
  b <- predict(m, gd, block = TRUE, se.fit = TRUE, local = cfg("maxmin"))$se.fit
  expect_identical(a, b)

  # grts is reproducible only under a fixed seed
  set.seed(11); g1 <- predict(m, gd, block = TRUE, se.fit = TRUE, local = cfg("grts"))$se.fit
  set.seed(11); g2 <- predict(m, gd, block = TRUE, se.fit = TRUE, local = cfg("grts"))$se.fit
  set.seed(12); g3 <- predict(m, gd, block = TRUE, se.fit = TRUE, local = cfg("grts"))$se.fit
  expect_identical(g1, g2)
  expect_false(isTRUE(all.equal(g1, g3)))

  # invalid ordering errors with the sprnorm()-style message
  expect_error(
    predict(m, gd, block = TRUE, se.fit = TRUE, local = cfg("not-an-ordering")),
    "ordering must be"
  )

  # unset ordering defaults to grts (seed-dependent), not maxmin (deterministic)
  cfg_default <- list(method = "all", method_new = "basis", size_new = 800)
  set.seed(21); d1 <- predict(m, gd, block = TRUE, se.fit = TRUE, local = cfg_default)$se.fit
  set.seed(22); d2 <- predict(m, gd, block = TRUE, se.fit = TRUE, local = cfg_default)$se.fit
  expect_false(isTRUE(all.equal(d1, d2)))

  # coincident coordinates still work under the default (grts handles them;
  # maxmin could not)
  gd_dup <- rbind(gd, gd[1:50, ])
  gd_dup$x <- sin(6 * gd_dup$xco) * cos(5 * gd_dup$yco)
  expect_error(
    predict(m, gd_dup, block = TRUE, se.fit = TRUE, local = cfg_default),
    NA
  )
})

test_that("block node ordering uses the anisotropy-transformed coordinates", {
  set.seed(7)
  n_obs <- 250
  d <- data.frame(xco = runif(n_obs), yco = runif(n_obs))
  d$s <- as.numeric(sprnorm(
    spcov_params("exponential", de = 0.9, ie = 0.1, range = 0.2, rotate = pi / 6, scale = 0.4),
    data = d, xcoord = xco, ycoord = yco
  ))
  d$z <- 2 + d$s
  m <- splm(z ~ 1, d,
    spcov_initial = spcov_initial("exponential", de = 0.9, ie = 0.1, range = 0.2, rotate = pi / 6, scale = 0.4, known = "given"),
    xcoord = xco, ycoord = yco, anisotropy = TRUE
  )
  gd <- bp_grid(4000)

  # exact path unaffected by the ordering; approximate path must still run be close
  ex <- predict(m, gd, block = TRUE, se.fit = TRUE, local = FALSE)
  set.seed(1)
  ap <- predict(m, gd, block = TRUE, se.fit = TRUE,
    local = list(method = "all", method_new = "basis", size_new = 1500, ordering = "grts"))
  expect_equal(ap$fit, ex$fit, tolerance = 1e-8)
  expect_equal(ap$se.fit, ex$se.fit, tolerance = 0.05)

  # Raw prediction coordinates must reach covmatrix(); transforming the grid
  # before this point would apply anisotropy a second time.
  gd_reference <- gd[seq_len(12), , drop = FALSE]
  S <- as.matrix(covmatrix(m))
  C <- as.matrix(covmatrix(m, gd_reference, cov_type = "pred.obs"))
  S0 <- as.matrix(covmatrix(m, gd_reference, cov_type = "pred.pred"))
  X <- model.matrix(m)
  X0 <- model.matrix(delete.response(terms(m)), gd_reference)
  V <- vcov(m)
  Q <- solve(S)
  B <- V %*% t(X) %*% Q
  W <- X0 %*% B + C %*% Q %*% (diag(nrow(X)) - X %*% B)
  G <- X0 - C %*% Q %*% X
  pred_cov <- S0 - C %*% Q %*% t(C) + G %*% V %*% t(G)
  expected_fit <- mean(W %*% model.response(model.frame(m)))
  expected_se <- sqrt(sum(pred_cov)) / nrow(gd_reference)
  got <- predict(m, gd_reference, block = TRUE, se.fit = TRUE, local = FALSE)
  expect_equal(got$fit, expected_fit, tolerance = 1e-8)
  expect_equal(got$se.fit, expected_se, tolerance = 1e-8)
})

test_that("subset block variance remains positive", {
  d <- data.frame(cx = seq(0, 1, length.out = 40), cy = 0, z = sin(seq_len(40)))
  fit <- splm(z ~ 1, d, xcoord = cx, ycoord = cy, spcov_initial =
    spcov_initial("exponential", de = .8, ie = .01, range = 1/3, known = "given"))
  nd <- d[seq(1, 40, by = 2), ]
  local <- list(method = "all", method_new = "subset", size_new = 4, ordering = "none")
  p <- predict(fit, nd, block = TRUE, se.fit = TRUE, local = local)
  expect_equal(p$se.fit^2, .002723610935126, tolerance = 1e-10)
  expect_gte(p$se.fit^2, .01 / NROW(nd))
  local$size_new <- NROW(nd)
  expect_identical(predict(fit, nd, block = TRUE, se.fit = TRUE, local = local),
    predict(fit, nd, block = TRUE, se.fit = TRUE, local = FALSE))
})

block_index_reference <- function(fit, nd, keep, nodes = seq_len(nrow(nd)), method_new = "subset") {
  cov_new <- covmatrix(fit, nd, cov_type = "pred.pred")
  cov_cross <- covmatrix(fit, nd, cov_type = "pred.obs")
  if (method_new == "basis") {
    c0 <- colMeans(cov_cross[, keep, drop = FALSE])
    s0 <- mean(cov_new[, nodes, drop = FALSE])
  } else {
    c0 <- colMeans(cov_cross[nodes, keep, drop = FALSE])
    nugget <- spcov_ie_stabilized(coef(fit, "spcov"), fit$diagtol)
    s0 <- mean(cov_new[nodes, nodes, drop = FALSE]) + nugget * (1/nrow(nd) - 1/length(nodes))
  }
  w <- solve(covmatrix(fit)[keep, keep, drop = FALSE], c0)
  X <- model.matrix(fit)[keep, , drop = FALSE]
  mat <- get_newdata_model_matrix(fit, nd)
  x0 <- colMeans(mat$newdata_model)
  h <- x0 - drop(crossprod(X, w))
  frame <- model.frame(fit)
  y <- model.response(frame)
  if (!is.null(model.offset(frame))) y <- y - model.offset(frame)
  offset <- if (is.null(mat$offset)) 0 else mean(mat$offset)
  list(fit = offset + sum(x0 * coef(fit)) + sum(w * (y[keep] - drop(X %*% coef(fit)))),
    se.fit = sqrt(max(0, s0 - sum(w*c0) + drop(crossprod(h, vcov(fit) %*% h)))))
}

test_that("block prediction uses the same rows for explicit subsets", {
  set.seed(73)
  d <- data.frame(cx = runif(36), cy = runif(36), x = rnorm(36), z = rnorm(36),
    group = factor(rep(letters[1:3], each = 12)), part = factor(rep(1:2, 18)))
  nd <- d[1:6, ]
  keep <- c(1L, 2L, 9L, 13L, 14L)
  for (structure in c("plain", "random", "partition")) {
    fit <- splm(z ~ x, d, xcoord = cx, ycoord = cy,
      random = if (structure == "random") ~group else NULL,
      randcov_initial = if (structure == "random") randcov_initial(group = .3, known = "given") else NULL,
      partition_factor = if (structure == "partition") ~part else NULL,
      spcov_initial = spcov_initial("exponential", de = .8, ie = .2, range = .3, known = "given"))
    subset_fit <- fit
    subset_fit$obdata <- fit$obdata[rev(keep), , drop = FALSE]
    expect_equal(unname(covmatrix(subset_fit)), unname(covmatrix(fit)[rev(keep), rev(keep)]), tolerance = 1e-12)
    ref <- block_index_reference(fit, nd, keep, c(2L, 5L))
    for (method in c("all", "covariance", "distance")) {
      actual <- predict(fit, nd, block = TRUE, se.fit = TRUE, local = list(
        method = method, size = 1, index = rev(keep), method_new = "subset",
        index_new = c(5, 2), size_new = Inf, ordering = "none"))
      expect_equal(actual, ref, tolerance = 1e-12)
    }
  }
})

test_that("block subset positions follow fitted and stored sf rows", {
  set.seed(74)
  d <- data.frame(cx = runif(30), cy = runif(30), x = rnorm(30), off = rnorm(30), z = rnorm(30))
  missing <- c(2L, 12L, 22L)
  d$z[missing] <- NA
  sf_data <- sf::st_as_sf(d, coords = c("cx", "cy"))
  fit <- splm(z ~ x + offset(off), sf_data, local = list(index = rep(c(3, 1, 2), 9)),
    spcov_initial = spcov_initial("exponential", de = .8, ie = .2, range = .3, known = "given"))
  nd <- sf_to_df(sf_data[missing, ])
  keep <- c(2L, 3L, 7L, 11L)
  expect_identical(rownames(model.frame(fit)), rownames(d)[!is.na(d$z)])
  ref <- block_index_reference(fit, nd, keep, 2L)
  local <- list(method = "all", index = keep, method_new = "subset", index_new = 2L)
  baseline <- predict(fit, block = TRUE, se.fit = TRUE, local = local)
  expect_equal(baseline, ref, tolerance = 1e-12)
  local$index <- seq_len(fit$n) %in% keep
  local$index_new <- c(FALSE, TRUE, FALSE)
  expect_identical(predict(fit, sf_data[missing, ], block = TRUE, se.fit = TRUE, local = local), baseline)
  local$index <- rev(keep)
  expect_identical(predict(fit, block = TRUE, se.fit = TRUE, local = local), baseline)
  expect_identical(predict(fit, nd, block = TRUE, type = "terms", local = local),
    predict(fit, nd, block = TRUE, type = "terms"))
  expect_identical(predict(fit, block = TRUE, interval = "confidence", local = local),
    predict(fit, block = TRUE, interval = "confidence"))
})

test_that("explicit grid nodes replace selection with the shipped method defaults", {
  fit <- bp_model(35, de = .8, ie = .2, range = .3)
  nd <- bp_grid(20)
  set.seed(80)
  nodes <- get_decorrelate_order("grts", nd$xco, nd$yco)$order[1:4]
  for (method_new in c("basis", "subset")) {
    set.seed(80)
    expected <- predict(fit, nd, block = TRUE, se.fit = TRUE,
      local = list(method_new = method_new, size_new = 4))
    actual <- predict(fit, nd, block = TRUE, se.fit = TRUE,
      local = list(method_new = method_new, index_new = nodes, ordering = "none"))
    expect_equal(actual, expected, tolerance = 1e-10)
  }
  expect_identical(predict(fit, nd, block = TRUE, se.fit = TRUE,
    local = list(method = "all", index_new = rev(seq_len(nrow(nd))), size_new = 1)),
    predict(fit, nd, block = TRUE, se.fit = TRUE, local = FALSE))
  for (name in c("index", "index_new")) {
    for (invalid in list(FALSE, c(1, NA), 0, 100, 1.5, Inf, c(1,1), integer(),
      "1", factor(1), matrix(1))) {
      expect_error(predict(fit, nd, block = TRUE, local = setNames(list(invalid), name)),
        paste0("local\\$", name))
    }
  }
  expect_error(predict(fit, nd, block = TRUE, local = list(index_new = 1, size_new = -1)), "size_new")
  expect_error(predict(fit, nd, block = TRUE, local = list(index_new = 1, ordering = "bad")), "ordering")
})

test_that("dense subset prediction retains the independent variance floor", {
  fit <- bp_model(2000, de = .8, ie = .01, range = .3)
  nd <- bp_grid(100)
  p <- predict(fit, nd, block = TRUE, se.fit = TRUE,
    local = list(method = "all", method_new = "subset", size_new = 25, ordering = "none"))
  expect_gte(p$se.fit^2, .01 / nrow(nd) * (1 - 1e-8))
})

test_that("subset covariance preserves anisotropy and the nugget floor", {
  set.seed(81)
  d <- data.frame(cx = runif(30), cy = runif(30), z = rnorm(30),
    group = factor(rep(letters[1:3], each = 10)), part = factor(rep(1:2, 15)))
  fit <- splm(z ~ 1, d, xcoord = cx, ycoord = cy, anisotropy = TRUE,
    partition_factor = ~part, random = ~group,
    randcov_initial = randcov_initial(group = .3, known = "given"),
    spcov_initial = spcov_initial("exponential", de = .8, ie = 0, range = .3,
      rotate = pi/6, scale = .4, known = "given"))
  fit$diagtol <- .012
  fit_1d <- splm(z ~ 1, d, xcoord = cx,
    spcov_initial = spcov_initial("exponential", de = .8, ie = .2, range = .3, known = "given"))
  fit_0d <- splm(z ~ 1, d, spcov_type = "ie")
  for (model in list(fit, fit_1d, fit_0d)) {
    for (keep in list(c(18L, 2L, 12L, 4L), 2L)) {
      selected <- model
      selected$obdata <- model$obdata[keep, , drop = FALSE]
      expect_equal(unname(covmatrix(selected)), unname(covmatrix(model)[keep, keep, drop = FALSE]),
        tolerance = 1e-12)
    }
  }
})
