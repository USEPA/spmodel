skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

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

test_that("method_new = 'subset' converges to the exact block prediction and variance", {
  gd <- bp_grid(6000, lims = c(0.1, 0.4)) # sub-region block
  m <- bp_model(220, de = 0.8, ie = 0.2, range = 0.15, formula = z ~ x)
  ex <- predict(m, gd, block = TRUE, se.fit = TRUE, local = FALSE)
  set.seed(4)
  sm <- predict(m, gd, block = TRUE, se.fit = TRUE, local = list(method = "all", method_new = "subset", size_new = 2000, ordering = "grts"))
  expect_equal(sm$fit, ex$fit, tolerance = 5e-3) # small block-mean error
  expect_equal(sm$se.fit, ex$se.fit, tolerance = 0.05)

  # the 1/size_new -> 1/G diagonal re-weight targets the exact (diagonal-counted) s0
  s0_today <- mean(covmatrix(m, newdata = gd, cov_type = "pred.pred"))
  de_ie <- sum(coef(m, "spcov")[c("de", "ie")])
  set.seed(4)
  nodes <- get_decorrelate_order("grts", gd$xco, gd$yco)$order[seq_len(2000)]
  raw <- spmodel:::get_block_quantities(m, gd[nodes, , drop = FALSE], nodes = seq_len(2000))$s0
  corrected <- raw + de_ie * (1 / nrow(gd) - 1 / 2000)
  expect_lt(abs(corrected - s0_today), 2e-3)
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
})
