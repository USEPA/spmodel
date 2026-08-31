skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

test_that("offsets are handled consistently for splm and spautor", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  # delta is deliberately not in the column space of cbind(1, x)
  exdata$delta <- sin(seq_len(nrow(exdata)))
  exdata$off <- log(2 + exdata$xcoord^2)
  exdata$cst <- rep(0.6, nrow(exdata))
  exdata$prop <- 0.6 * exdata$x
  exdata$yshift <- exdata$y + exdata$delta
  exdata$offshift <- exdata$off + exdata$delta
  newexdata$delta <- sin(seq_len(nrow(newexdata)) + 0.5)
  newexdata$off <- log(2 + newexdata$xcoord^2)
  newexdata$cst <- rep(0.6, nrow(newexdata))
  newexdata$prop <- 0.6 * newexdata$x
  newexdata$offshift <- newexdata$off + newexdata$delta

  fit <- function(form) {
    splm(form, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  }
  m_none <- fit(y ~ x)
  m_cst <- fit(y ~ x + offset(cst))
  m_prop <- fit(y ~ x + offset(prop))

  # (I) and (II): only the named coefficient moves
  expect_equal(unname(coef(m_cst) - coef(m_none)), c(-0.6, 0))
  expect_equal(unname(coef(m_prop) - coef(m_none)), c(0, -0.6))
  for (m in list(m_cst, m_prop)) {
    expect_equal(coef(m, type = "spcov"), coef(m_none, type = "spcov"))
    expect_equal(as.numeric(logLik(m)), as.numeric(logLik(m_none)))
    expect_equal(fitted(m), fitted(m_none))
    expect_equal(residuals(m), residuals(m_none))
    expect_equal(hatvalues(m), hatvalues(m_none))
    expect_equal(AIC(m), AIC(m_none))
    # every prediction and cross-validation path
    expect_equal(predict(m, newexdata), predict(m_none, newexdata))
    expect_equal(
      predict(m, newexdata, se.fit = TRUE)$se.fit,
      predict(m_none, newexdata, se.fit = TRUE)$se.fit
    )
    expect_equal(loocv(m, cv_predict = TRUE)$cv_predict, loocv(m_none, cv_predict = TRUE)$cv_predict)
    fi <- rep(1:4, length.out = nrow(exdata))
    expect_equal(
      kcv(m, folds_index = fi, cv_predict = TRUE)$cv_predict,
      kcv(m_none, folds_index = fi, cv_predict = TRUE)$cv_predict
    )
  }

  # (III) with a general offset: shifting y and the offset together must leave
  # every estimate alone and shift every prediction by exactly delta. An offset
  # left in the quantity being kriged would instead be partly absorbed by
  # betahat, and this is the check that sees it.
  a <- fit(y ~ x + offset(off))
  b <- fit(yshift ~ x + offset(offshift))
  expect_equal(coef(b), coef(a))
  expect_equal(coef(b, type = "spcov"), coef(a, type = "spcov"))
  expect_equal(residuals(b), residuals(a))
  expect_equal(
    as.numeric(predict(b, newexdata) - predict(a, newexdata)),
    newexdata$delta
  )
  expect_equal(
    as.numeric(loocv(b, cv_predict = TRUE)$cv_predict - loocv(a, cv_predict = TRUE)$cv_predict),
    exdata$delta
  )
  fi <- rep(1:4, length.out = nrow(exdata))
  expect_equal(
    as.numeric(kcv(b, folds_index = fi, cv_predict = TRUE)$cv_predict -
      kcv(a, folds_index = fi, cv_predict = TRUE)$cv_predict),
    exdata$delta
  )
  for (loc in list(FALSE, list(approximation = "vecchia"))) {
    set.seed(2)
    ca <- conditional(a, newexdata, samples = 5, local = loc)
    set.seed(2)
    cb <- conditional(b, newexdata, samples = 5, local = loc)
    expect_equal(as.numeric(cb - ca), rep(newexdata$delta, times = 5))
  }

  # spautor: prediction locations are the NA rows of the same data
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  np <- nrow(exdata_poly)
  exdata_poly$delta <- sin(seq_len(np))
  exdata_poly$off <- log(2 + exdata_poly$x^2)
  exdata_poly$offshift <- exdata_poly$off + exdata_poly$delta
  exdata_poly$yshift <- exdata_poly$y + exdata_poly$delta
  exdata_poly$y[1:5] <- NA
  exdata_poly$yshift[1:5] <- NA

  pa <- spautor(y ~ x + offset(off), exdata_poly, spcov_type = "car")
  pb <- spautor(yshift ~ x + offset(offshift), exdata_poly, spcov_type = "car")
  expect_equal(coef(pb), coef(pa))
  expect_equal(as.numeric(predict(pb) - predict(pa)), exdata_poly$delta[1:5])
  expect_equal(predict(pb, se.fit = TRUE)$se.fit, predict(pa, se.fit = TRUE)$se.fit)
})

test_that("offsets are handled consistently for spglm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  n <- nrow(exdata)
  exdata$count <- as.integer(round(abs(exdata$y) * 2))
  exdata$cst <- rep(0.6, n)
  exdata$prop <- 0.6 * exdata$x
  newexdata$cst <- rep(0.6, nrow(newexdata))
  newexdata$prop <- 0.6 * newexdata$x

  fit <- function(form, local = FALSE) {
    spglm(form, family = "poisson", exdata,
      spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = local
    )
  }
  m_none <- fit(count ~ x)
  m_cst <- fit(count ~ x + offset(cst))
  m_prop <- fit(count ~ x + offset(prop))

  expect_equal(unname(coef(m_cst) - coef(m_none)), c(-0.6, 0), tolerance = 1e-5)
  expect_equal(unname(coef(m_prop) - coef(m_none)), c(0, -0.6), tolerance = 1e-5)

  # m_prop is the informative one for anything evaluated at the fitted mean: the
  # offset varies across observations, so a leverage weight or residual computed
  # from the offset-free latent vector rather than the linear predictor differs.
  for (m in list(m_cst, m_prop)) {
    expect_equal(coef(m, type = "spcov"), coef(m_none, type = "spcov"), tolerance = 1e-5)
    expect_equal(as.numeric(logLik(m)), as.numeric(logLik(m_none)), tolerance = 1e-3)
    expect_equal(fitted(m), fitted(m_none), tolerance = 1e-5)
    expect_equal(hatvalues(m), hatvalues(m_none), tolerance = 1e-5)
    for (rt in c("response", "deviance", "pearson", "standardized")) {
      expect_equal(residuals(m, type = rt), residuals(m_none, type = rt), tolerance = 1e-5)
    }
    expect_equal(cooks.distance(m), cooks.distance(m_none), tolerance = 1e-5)
    expect_equal(predict(m, newexdata), predict(m_none, newexdata), tolerance = 1e-5)
    expect_equal(
      predict(m, newexdata, se.fit = TRUE)$se.fit,
      predict(m_none, newexdata, se.fit = TRUE)$se.fit,
      tolerance = 1e-5
    )
    # the nearest-neighbour prediction path evaluates the latent-process
    # variance adjustment separately from the dense path
    loc <- list(method = "covariance", size = 30)
    expect_equal(
      predict(m, newexdata, se.fit = TRUE, local = loc)$se.fit,
      predict(m_none, newexdata, se.fit = TRUE, local = loc)$se.fit,
      tolerance = 1e-5
    )
    expect_equal(
      loocv(m, cv_predict = TRUE)$cv_predict,
      loocv(m_none, cv_predict = TRUE)$cv_predict,
      tolerance = 1e-4
    )
    fi <- rep(1:4, length.out = n)
    expect_equal(
      kcv(m, folds_index = fi, cv_predict = TRUE)$cv_predict,
      kcv(m_none, folds_index = fi, cv_predict = TRUE)$cv_predict,
      tolerance = 1e-4
    )
    set.seed(1)
    c_m <- conditional(m, newexdata, samples = 5)
    set.seed(1)
    c_none <- conditional(m_none, newexdata, samples = 5)
    expect_equal(c_m, c_none, tolerance = 1e-4)
    # the vecchia approximation reaches the latent-process variance adjustment
    # by a different route than the default low-rank one, so it is checked too
    vec <- list(approximation = "vecchia")
    set.seed(2)
    v_m <- conditional(m, newexdata, samples = 5, local = vec)
    set.seed(2)
    v_none <- conditional(m_none, newexdata, samples = 5, local = vec)
    expect_equal(v_m, v_none, tolerance = 1e-4)
  }

  # local (big data) fitting partitions the data, and the offset has to survive
  # into every partition's contribution to the latent-vector solve
  loc <- list(index = rep(1:3, length.out = n), var_adjust = "none")
  l_none <- fit(count ~ x, local = loc)
  l_cst <- fit(count ~ x + offset(cst), local = loc)
  l_prop <- fit(count ~ x + offset(prop), local = loc)
  expect_equal(unname(coef(l_cst) - coef(l_none)), c(-0.6, 0), tolerance = 1e-5)
  expect_equal(unname(coef(l_prop) - coef(l_none)), c(0, -0.6), tolerance = 1e-5)
  expect_equal(fitted(l_cst), fitted(l_none), tolerance = 1e-5)
  expect_equal(fitted(l_prop), fitted(l_none), tolerance = 1e-5)

  # A general offset -- one that is neither constant nor a linear combination of
  # the columns of X, so betahat cannot absorb it and invariances (I) and (II)
  # cannot see it. Checked against universal kriging written out directly: the
  # offset comes out of the latent vector, the remaining n - 1 observations
  # predict the held-out one, and that row's own offset goes back on.
  exdata$off <- log(2 + exdata$xcoord^2)
  m_off <- fit(count ~ x + offset(off))
  Sig <- as.matrix(covmatrix(m_off))
  X <- model.matrix(m_off)
  w <- as.numeric(fitted(m_off, type = "link")) - exdata$off
  cv_ref <- vapply(seq_len(n), function(i) {
    Si <- chol2inv(chol(Sig[-i, -i]))
    Xm <- X[-i, , drop = FALSE]
    wm <- w[-i]
    b <- solve(crossprod(Xm, Si %*% Xm), crossprod(Xm, Si %*% wm))
    as.numeric(X[i, , drop = FALSE] %*% b +
      crossprod(Sig[-i, i], Si %*% (wm - Xm %*% b))) + exdata$off[i]
  }, numeric(1))
  expect_equal(loocv(m_off, cv_predict = TRUE, type = "link")$cv_predict, cv_ref)
})

test_that("offsets are handled consistently for spgautor", {
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

  np <- nrow(exdata_poly)
  exdata_poly$count <- as.integer(round(abs(exdata_poly$y) * 2))
  exdata_poly$cst <- rep(0.6, np)
  exdata_poly$prop <- 0.6 * exdata_poly$x
  exdata_poly$count[1:5] <- NA

  fit <- function(form) spgautor(form, family = "poisson", exdata_poly, spcov_type = "car")
  m_none <- fit(count ~ x)
  m_cst <- fit(count ~ x + offset(cst))
  m_prop <- fit(count ~ x + offset(prop))

  expect_equal(unname(coef(m_cst) - coef(m_none)), c(-0.6, 0), tolerance = 1e-5)
  expect_equal(unname(coef(m_prop) - coef(m_none)), c(0, -0.6), tolerance = 1e-5)
  for (m in list(m_cst, m_prop)) {
    expect_equal(fitted(m), fitted(m_none), tolerance = 1e-5)
    expect_equal(hatvalues(m), hatvalues(m_none), tolerance = 1e-5)
    expect_equal(residuals(m, type = "standardized"), residuals(m_none, type = "standardized"), tolerance = 1e-5)
    expect_equal(predict(m), predict(m_none), tolerance = 1e-5)
    expect_equal(predict(m, se.fit = TRUE)$se.fit, predict(m_none, se.fit = TRUE)$se.fit, tolerance = 1e-5)
    expect_equal(
      loocv(m, cv_predict = TRUE)$cv_predict,
      loocv(m_none, cv_predict = TRUE)$cv_predict,
      tolerance = 1e-4
    )
  }
})

test_that("offsets are handled consistently in block prediction", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  exdata$cst <- rep(0.6, nrow(exdata))
  exdata$off <- log(2 + exdata$xcoord^2)
  newexdata$cst <- rep(0.6, nrow(newexdata))
  newexdata$off <- log(2 + newexdata$xcoord^2)

  fit <- function(form) {
    splm(form, exdata, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  }
  m_none <- fit(y ~ x)
  m_cst <- fit(y ~ x + offset(cst))
  m_off <- fit(y ~ x + offset(off))

  # a block prediction is a single value regardless of the offset -- the
  # observed-data offset must not leak into the returned vector
  b_none <- predict(m_none, newexdata, block = TRUE, se.fit = TRUE)
  b_cst <- predict(m_cst, newexdata, block = TRUE, se.fit = TRUE)
  expect_length(b_cst$fit, 1)
  expect_equal(b_cst$fit, b_none$fit)
  expect_equal(b_cst$se.fit, b_none$se.fit)

  # the block's own offset is the average of the newdata offsets, because a
  # block prediction is the average of the point predictions
  b_off <- predict(m_off, newexdata, block = TRUE)
  expect_length(b_off, 1)
  expect_equal(as.numeric(b_off), mean(predict(m_off, newexdata)))

  expect_equal(nrow(predict(m_off, newexdata, block = TRUE, interval = "confidence")), 1)
  expect_equal(nrow(predict(m_off, newexdata, block = TRUE, interval = "prediction")), 1)
})

# An offset is an explanatory variable whose coefficient is fixed at one rather
# than estimated, so it is a known shift that must never be incorporated into betahat.
# Three exact algebraic invariances follow, and between them they pin down every
# place in the package that has to handle an offset:
#
#  (I)   offset = a constant c, in a model with an intercept: the intercept
#        shifts by exactly -c and nothing else changes at all.
#  (II)  offset = c * x for an x already in the model: the coefficient on x
#        shifts by exactly -c and nothing else changes. Unlike (I) the offset
#        now varies across observations, which is what makes it able to see
#        quantities evaluated at the fitted mean (leverage weights, residuals).
#  (III) (a not very relevant case) shifting both the response and the offset by the same arbitrary vector
#        delta leaves the offset-free problem (y + delta) - (o + delta) = y - o
#        untouched, so every estimate is unchanged and every prediction shifts
#        by exactly delta. This is the only one of the three that admits a delta
#        outside the column space of X, and so the only one that can catch code
#        which lets the fixed effects incorporate part of the offset.