skip_on_cran()
skip_if_not(
  identical(Sys.getenv("SPMODEL_RUN_EXTRAS"), "true"),
  "set Sys.setenv(SPMODEL_RUN_EXTRAS = 'true') before devtools::test() to run the extras suite"
)

set.seed(1)

# SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_M.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

# add variables
n <- NROW(exdata)
exdata$bern <- rbinom(n, size = 1, prob = 0.5)
exdata$bernfac <- factor(ifelse(exdata$bern == 0, "a", "b"))
exdata$size <- 5
exdata$bin <- rbinom(n, size = exdata$size, prob = 0.5)
exdata$prop <- runif(n)
exdata$count <- rpois(n, lambda = 5)
exdata$cont <- rgamma(n, shape = 1, rate = 1)
exdata$offset <- 1.2

spglm_prediction_fixture <- function() {
  set.seed(90210)
  n <- 30
  dat <- data.frame(
    cx = rep(seq_len(6), 5),
    cy = rep(seq_len(5), each = 6),
    x = seq(-1, 1, length.out = n),
    off = 0.2 * sin(seq_len(n)),
    part = factor(rep(c("a", "b", "c"), length.out = n))
  )
  dat$pois <- rpois(n, exp(-0.2 + 0.3 * dat$x + dat$off))
  pred <- data.frame(
    cx = c(1.2, 2.4, 3.6, 4.8, 2.1, 5.2),
    cy = c(1.4, 2.7, 3.3, 4.2, 4.5, 1.8),
    x = seq(-0.8, 0.8, length.out = 6),
    off = 0.1 * cos(seq_len(6)),
    part = factor(c("a", "b", "c", "a", "b", "c"), levels = c("a", "b", "c"))
  )
  list(data = dat, pred = pred)
}

spglm_prediction_reference <- function(object, newdata) {
  S <- as.matrix(covmatrix(object))
  C <- as.matrix(covmatrix(object, newdata, cov_type = "pred.obs"))
  S0 <- as.matrix(covmatrix(object, newdata, cov_type = "pred.pred"))
  X <- model.matrix(object)
  X0 <- model.matrix(delete.response(terms(object)), newdata)
  V <- vcov(object, var_correct = FALSE)
  Q <- solve(S)
  B <- V %*% t(X) %*% Q
  W <- X0 %*% B + C %*% Q %*% (diag(nrow(X)) - X %*% B)
  G <- X0 - C %*% Q %*% X
  variance <- S0 - C %*% Q %*% t(C) + G %*% V %*% t(G)
  P <- Q - Q %*% X %*% B
  variance <- variance + W %*% solve(P + diag(exp(fitted(object, "link")))) %*% t(W)
  w <- fitted(object, "link") - model.offset(model.frame(object))
  list(fit = as.vector(newdata$off + W %*% w), variance = variance)
}

test_that("the model runs for binomial data", {
  spgmod <- spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  expect_s3_class(spgmod, "spglm")
  expect_vector(AUROC(spgmod))
  spgmod <- spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  expect_s3_class(spgmod, "spglm")
  expect_vector(AUROC(spgmod))
  expect_error(spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
  spgmod <- suppressWarnings(spglm(cbind(bin, size) ~ x, family = "binomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "ml")) # ml boundary warning
  expect_s3_class(spgmod, "spglm")
  expect_error(AUROC(spgmod))
  spgmod <- spglm(y > 0 ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  expect_s3_class(spgmod, "spglm")
  expect_vector(AUROC(spgmod))


  # complicated models
  # optim non-convergence warning expected for these numerically marginal fixtures
  expect_error(suppressWarnings(spglm(bern ~ x + offset(offset), family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE)), NA)
  # check when data have an unequal number of local observations
  expect_error(spglm(bern ~ x + offset(offset), family = binomial, data = exdata[-1, , drop = FALSE], xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE), NA)
  expect_error(suppressWarnings(spglm(bernfac ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, anisotropy = TRUE, local = TRUE)), NA)
  expect_error(suppressWarnings(spglm(cbind(bin, size) ~ x, family = "binomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "ml", partition_factor = ~group)), NA) # ml boundary warning

  # check glm
  spgmod1 <- spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none")
  spgmod2 <- spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential")
  expect_false(isTRUE(all.equal(as.vector(coef(spgmod1)), as.vector(coef(spgmod2)), tolerance = 0.0001)))
  gmod1 <- stats::glm(bern ~ x, family = binomial, data = exdata)
  expect_true(isTRUE(all.equal(as.vector(coef(spgmod1)), as.vector(coef(gmod1)), tolerance = 0.0001)))
})

test_that("the model runs for proportion data", {
  expect_error(spglm(prop ~ x, family = "beta", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
  expect_error(spglm(prop ~ x, family = beta, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"), NA)

  # complicated models
  expect_error(spglm(prop ~ x, family = "beta", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", anisotropy = TRUE), NA)
  expect_error(spglm(prop ~ x + offset(offset),
    family = beta, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "ie", estmethod = "reml",
    local = list(method = "kmeans")
  ), NA)
})

test_that("the model runs for count data", {
  spgmod <- spglm(count ~ x, family = poisson, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml")
  expect_s3_class(spgmod, "spglm")
  expect_error(AUROC(spgmod))
  expect_error(suppressWarnings(spglm(count ~ x, family = "nbinomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml")), NA) # ml boundary warning

  # complicated models
  expect_error(spglm(count ~ x + offset(offset),
    family = poisson, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml",
    partition_factor = ~group, local = TRUE
  ), NA)
  expect_error(suppressWarnings(spglm(count ~ x, family = "nbinomial", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml", random = ~group)), NA) # ml boundary warning
})

test_that("the model runs for continuous data", {
  expect_error(spglm(cont ~ x, family = "Gamma", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml"), NA)
  expect_error(spglm(cont ~ x, family = inverse.gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"), NA)

  # SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
  expect_error(spglm(cont ~ x, family = gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"))


  # complicated models
  expect_error(spglm(cont ~ x + offset(offset),
    family = "Gamma", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "spherical", estmethod = "reml",
    random = ~ group + subgroup, local = TRUE
  ), NA)
  expect_error(spglm(cont ~ x, family = inverse.gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml", partition_factor = ~subgroup, anisotropy = TRUE), NA)
  ## SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
  expect_error(spglm(cont ~ x,
    family = gaussian, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml",
    random = ~group
  ))
})

test_that("the model runs on other data sets", {
  expect_error(spglm(abs(y) ~ x, family = "Gamma", data = exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"), NA)
  expect_error(suppressWarnings(spglm(abs(y) ~ x, family = Gamma, data = exdata_poly, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml")), NA) # POINT, xcoord, ycoord warnings


  # complicated models
  expect_error(spglm(abs(y) ~ x, family = "Gamma", data = exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group, partition_factor = ~group, anisotropy = TRUE, local = TRUE), NA)
  expect_error(suppressWarnings(spglm(abs(y) ~ x,
    family = Gamma, data = exdata_poly, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml",
    random = ~group, partition_factor = ~subgroup
  )), NA) # POINT, xcoord, ycoord warnings
})


test_that("emmeans works", {
  spcov_type <- "exponential"
  spgmod <- spglm(abs(y) ~ x * group, family = "Gamma", exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(as.matrix(model.frame(delete.response(terms(spgmod)), spgmod$obdata)), as.matrix(emmeans::recover_data(spgmod)))
  expect_error(emmeans::emmeans(spgmod, ~ group, by = "x"), NA)
})

test_that("emmeans works missing", {
  spcov_type <- "exponential"
  spgmod <- spglm(abs(y) ~ x * group, family = "Gamma", exdata_M, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml")
  expect_equal(as.matrix(model.frame(delete.response(terms(spgmod)), spgmod$obdata)), as.matrix(emmeans::recover_data(spgmod)))
  expect_error(emmeans::emmeans(spgmod, ~ group, by = "x"), NA)
})

test_that("range_constrain works", {
  spcov_type <- "exponential"
  expect_error(spglm(abs(y) ~ x * group, family = "Gamma", exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "reml", range_constrain = TRUE), NA)
  expect_error(suppressWarnings(spglm(abs(y) ~ x * group, family = "Gamma", exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = spcov_type, estmethod = "ml", range_constrain = TRUE)), NA) # ml boundary warning
})

test_that("ml boundary regression: an optimizer that wanders to the floor is reconciled and flagged", {
  # reproduces the scenario documented in the ml-boundary write-up: data simulated
  # with genuine spatial structure, whose exponential ml fit still converges to a
  # tiny ie chasing the boundary artifact
  set.seed(42)
  n_side <- 12
  coords <- expand.grid(xcoord = seq_len(n_side), ycoord = seq_len(n_side))
  coords$x <- rnorm(nrow(coords))

  sp_params <- spcov_params("exponential", de = 0.6, ie = 0.1, range = 4, rotate = 0, scale = 1)
  coords$y <- sprpois(sp_params, mean = 0.5 + 0.3 * coords$x, data = coords, xcoord = xcoord, ycoord = ycoord)

  expect_warning(
    mod_exp <- spglm(y ~ x, family = "poisson", data = coords, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"),
    NA # de stays substantial (see the write-up), so this specific fit should not trip the boundary warning
  )
  # the optimizer wanders below the floor for ie; the reported value must be reconciled to it
  expect_true(coef(mod_exp, type = "spcov")[["ie"]] >= 1e-4)

  expect_warning(
    mod_none <- spglm(y ~ x, family = "poisson", data = coords, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "ml"),
    "numerical boundary"
  )
  # reml should not show the same AIC preference for "none" over "exponential"
  mod_exp_reml <- spglm(y ~ x, family = "poisson", data = coords, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  mod_none_reml <- spglm(y ~ x, family = "poisson", data = coords, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml")
  expect_true(AIC(mod_exp_reml) < AIC(mod_none_reml))
})

test_that("anisotropic prediction uncertainty uses raw coordinates", {
  fixture <- spglm_prediction_fixture()
  dat <- fixture$data
  pred <- fixture$pred
  fit <- spglm(pois ~ x + cx + offset(off), "poisson", dat,
    xcoord = cx, ycoord = cy, anisotropy = TRUE,
    spcov_initial = spcov_initial("exponential",
      de = 0.15, ie = 0.12, range = 4, rotate = 0.6, scale = 0.45,
      known = "given"
    )
  )
  expected <- spglm_prediction_reference(fit, pred)
  got <- predict(fit, pred, type = "link", se.fit = TRUE, local = FALSE)
  expect_equal(unname(got$fit), expected$fit, tolerance = 1e-8)
  expect_equal(unname(got$se.fit), unname(sqrt(diag(expected$variance))), tolerance = 1e-8)
})

test_that("partitioned local prediction keeps row-dependent inputs aligned", {
  fixture <- spglm_prediction_fixture()
  dat <- fixture$data
  pred <- fixture$pred
  fit <- spglm(pois ~ x + offset(off), "poisson", dat,
    xcoord = cx, ycoord = cy, partition_factor = ~part,
    spcov_initial = spcov_initial("exponential",
      de = 0.15, ie = 0.12, range = 4, known = "given"
    )
  )
  exact <- predict(fit, pred, type = "link", se.fit = TRUE, local = FALSE)
  for (method in c("covariance", "distance")) {
    for (byrow_threshold in c(0, Inf)) {
      local <- predict(fit, pred, type = "link", se.fit = TRUE,
        local = list(
          method = method, size = nrow(dat), parallel = FALSE,
          byrow_threshold = byrow_threshold
        )
      )
      expect_equal(local$fit, exact$fit, tolerance = 1e-8)
      expect_equal(local$se.fit, exact$se.fit, tolerance = 1e-8)
    }
  }

  weights <- predict(fit, pred, type = "weight",
    local = list(method = "covariance", size = 5, parallel = FALSE)
  )
  expect_equal(ncol(weights), nrow(dat))
})

test_that("covariance components retain original observation order", {
  fixture <- spglm_prediction_fixture()
  dat <- fixture$data
  init <- spcov_initial("exponential", de = 0.15, ie = 0.12, range = 4, known = "given")
  fits <- list(
    partition = spglm(pois ~ x + offset(off), "poisson", dat,
      xcoord = cx, ycoord = cy, partition_factor = ~part, spcov_initial = init
    ),
    local = spglm(pois ~ x + offset(off), "poisson", dat,
      xcoord = cx, ycoord = cy, spcov_initial = init,
      local = list(index = rep(1:3, length.out = nrow(dat)), var_adjust = "theoretical")
    )
  )

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    r <- fitted(fit, "link") - model.offset(model.frame(fit)) -
      as.vector(model.matrix(fit) %*% coef(fit))
    S <- as.matrix(covmatrix(fit))
    if (nm == "local") {
      S <- S * outer(fit$local_index, fit$local_index, "==")
    }
    nugget <- coef(fit, "spcov")[["ie"]] * as.vector(solve(S, r))
    expect_equal(unname(fitted(fit, "spcov")$ie), nugget, tolerance = 1e-8)
    expect_named(fitted(fit, "spcov")$ie, as.character(fit$observed_index))
  }
})

test_that("spglm numerical floor is excluded from BLUPs and included in prediction variance", {
  n_obs <- 18
  dat <- data.frame(
    xc = seq(0, 17),
    yc = rep(c(0, 0.5), 9),
    x = seq(-1, 1, length.out = n_obs),
    off = 0.15 * sin(seq_len(n_obs)),
    y = c(1, 2, 1, 3, 2, 4, 2, 3, 5, 4, 6, 5, 7, 5, 8, 7, 9, 8)
  )
  newdata <- data.frame(
    xc = c(2.5, 8.5, 15.5),
    yc = c(0.25, 0.75, 0.25),
    x = c(-0.6, 0.1, 0.8),
    off = c(0.1, -0.05, 0.2)
  )
  de <- 0.4
  ie <- 0
  range <- 4
  fit <- spglm(y ~ x + offset(off), "poisson", dat,
    xcoord = xc, ycoord = yc,
    spcov_initial = spcov_initial("exponential",
      de = de, ie = ie, range = range, known = "given"
    )
  )

  dist_obs <- sqrt(
    outer(dat$xc, dat$xc, "-")^2 + outer(dat$yc, dat$yc, "-")^2
  )
  dist_pred <- sqrt(
    outer(newdata$xc, dat$xc, "-")^2 + outer(newdata$yc, dat$yc, "-")^2
  )
  K <- de * exp(-dist_obs / range)
  C <- de * exp(-dist_pred / range)
  effective_ie <- max(ie, 1e-4 * de, fit$diagtol)
  V <- K + diag(effective_ie, n_obs)
  V_inv <- solve(V)
  X <- model.matrix(fit)
  X0 <- model.matrix(delete.response(terms(fit)), newdata)
  beta <- coef(fit)
  w <- fitted(fit, "link") - model.offset(model.frame(fit))
  residual_weight <- V_inv %*% (w - X %*% beta)
  B <- vcov(fit, var_correct = FALSE)

  expect_equal(coef(fit, "spcov")[["ie"]], ie)
  expect_equal(
    unname(fitted(fit, "spcov")$de),
    as.numeric(K %*% residual_weight),
    tolerance = 1e-10
  )
  expect_equal(unname(fitted(fit, "spcov")$ie), rep(0, n_obs))

  H <- X0 - C %*% V_inv %*% X
  fit_reference <- as.numeric(X0 %*% beta + C %*% residual_weight + newdata$off)
  var_reference <- de + effective_ie - rowSums((C %*% V_inv) * C) +
    rowSums((H %*% B) * H)
  full_prediction <- predict(fit, newdata,
    type = "link", se.fit = TRUE, var_correct = FALSE, local = FALSE
  )
  expect_equal(unname(full_prediction$fit), fit_reference, tolerance = 1e-10)
  expect_equal(unname(full_prediction$se.fit), sqrt(as.numeric(var_reference)), tolerance = 1e-10)

  local_prediction <- predict(fit, newdata,
    type = "link", se.fit = TRUE, var_correct = FALSE,
    local = list(method = "distance", size = 8, parallel = FALSE)
  )
  local_reference <- lapply(seq_len(nrow(newdata)), function(i) {
    keep <- order(dist_pred[i, ])[seq_len(8)]
    V_local <- V[keep, keep, drop = FALSE]
    C_local <- C[i, keep, drop = FALSE]
    X_local <- X[keep, , drop = FALSE]
    w_local <- w[keep]
    V_local_inv <- solve(V_local)
    H_local <- X0[i, , drop = FALSE] - C_local %*% V_local_inv %*% X_local
    list(
      fit = as.numeric(
        X0[i, , drop = FALSE] %*% beta +
          C_local %*% V_local_inv %*% (w_local - X_local %*% beta) + newdata$off[i]
      ),
      var = as.numeric(
        de + effective_ie - C_local %*% V_local_inv %*% t(C_local) +
          H_local %*% B %*% t(H_local)
      )
    )
  })
  expect_equal(
    unname(local_prediction$fit),
    vapply(local_reference, `[[`, numeric(1), "fit"),
    tolerance = 1e-10
  )
  expect_equal(
    unname(local_prediction$se.fit),
    sqrt(vapply(local_reference, `[[`, numeric(1), "var")),
    tolerance = 1e-10
  )
})

test_that("size > 1 binomial spglm fit checks fitted probabilities, not fitted successes", {
  dat <- data.frame(
    x = seq(-1, 1, length.out = 30),
    size = rep(c(10, 20, 30), 10),
    xcoord = rep(1:6, 5),
    ycoord = rep(1:5, each = 6),
    off = 0.3 * sin(1:30)
  )
  dat$successes <- round(dat$size * (0.4 + 0.05 * dat$x))
  dat$failures <- dat$size - dat$successes
  dat$successes[c(3, 18)] <- NA

  expect_warning(
    fit <- spglm(cbind(successes, failures) ~ x + offset(off),
      family = "binomial", data = dat, xcoord = xcoord, ycoord = ycoord,
      spcov_initial = spcov_initial("exponential", de = 0.1, ie = 0.1, range = 1,
        known = c("de", "ie", "range")
      )
    ),
    NA
  )

  counts <- fitted(fit, type = "response")
  probabilities <- expit(fitted(fit, type = "link"))
  expect_true(all(counts > 1))
  expect_true(all(probabilities > 0.2 & probabilities < 0.6))
  expect_equal(length(probabilities), sum(!is.na(dat$successes)))
  expect_equal(counts, dat$size[!is.na(dat$successes)] * probabilities)
})

test_that("the binomial saturation warning retains its probability thresholds", {
  expect_warning(warn_fitted_saturation(rep(0.4, 100), "binomial"), NA)
  expect_warning(
    warn_fitted_saturation(c(rep(0, 49), rep(1, 49), 0.4, 0.6), "binomial"),
    NA
  )
  expect_warning(
    warn_fitted_saturation(c(rep(1e-7, 50), rep(1 - 1e-7, 49), 0.4), "binomial"),
    "Perfect separation"
  )
  expect_warning(warn_fitted_saturation(rep(8, 100), "poisson"), NA)
})

test_that("fitted-probability saturation warning fires under spatial separation but not on a sane binomial fit", {
  skip_if_not_installed("sf")

  # reproduces the scenario documented in the separation write-up: default
  # (reml) estmethod, de runs away, fitted probabilities fully saturate
  data(lake)
  lake_df <- sf::st_drop_geometry(lake)
  xy <- sf::st_coordinates(lake)
  lake_df$xcoord <- xy[, 1]
  lake_df$ycoord <- xy[, 2]
  lake_df$y <- as.numeric(lake_df$log_cond > 3)

  # this fixture also fails to converge (the saturation and non-convergence
  # diagnostics commonly co-occur); nest so both expected warnings are caught
  expect_warning(
    expect_warning(
      spglm(y ~ temp + elev, data = lake_df, family = binomial, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential"),
      "did not converge"
    ),
    "separation"
  )

  # a sane binomial fit (well below the saturation threshold, see the
  # calibration in the separation write-up) must not trip the warning
  expect_warning(
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml"),
    NA
  )
})

test_that("ml boundary warning fires for spcov_type = none but not reml", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "ml"),
    "numerical boundary"
  )
  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml"),
    NA
  )
})

test_that("ml boundary warning uses the actual fitted value and does not override a known ie", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  spinit <- spcov_initial("exponential", de = 1e-8, ie = 1e-8, range = 1, known = c("de", "ie", "range"))
  expect_warning(
    spmod <- spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_initial = spinit, estmethod = "ml"),
    "numerical boundary"
  )
  # a known (fixed) ie must never be reconciled to the diagtol floor
  expect_equal(coef(spmod, type = "spcov")[["ie"]], 1e-8)
})

test_that("no ml boundary warning for a well-identified spatial covariance", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "ml"),
    NA
  )
})

test_that("optim non-convergence warning fires for spglm", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", control = list(maxit = 1)),
    "did not converge"
  )
  expect_warning(
    spglm(y_pois ~ x, family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential"),
    NA
  )
})

test_that("anisotropy = TRUE with fully known covariance and dispersion parameters does not error", {
  # regression test for a bug where run_laploglik_dispatch_spglm() called
  # use_laploglik_known_anis() with an unnamed positional argument list one
  # short of its formals, silently dropping randcov_initial and erroring with
  # "argument \"randcov_initial\" is missing, with no default" -- for both the
  # no-random-effects case and the with-random-effects case
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  exdata$y_pois <- round(abs(exdata$y) * 3)

  spcov_initial_known <- spcov_initial("exponential",
                                       de = 1, ie = 0.5, range = 2, rotate = 0.3, scale = 0.7,
                                       known = "given"
  )
  dispersion_initial_known <- dispersion_initial("poisson", dispersion = 1, known = "given")

  # no random effects
  expect_error(
    spmod_no_randcov <- spglm(y_pois ~ x,
                              family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord,
                              spcov_initial = spcov_initial_known,
                              dispersion_initial = dispersion_initial_known,
                              anisotropy = TRUE
    ),
    NA
  )
  expect_equal(as.numeric(coef(spmod_no_randcov, type = "spcov")), as.numeric(spcov_initial_known$initial))

  # with a known random effect
  randcov_initial_known <- randcov_initial(group = 0.2, known = "given")
  expect_error(
    spmod_with_randcov <- spglm(y_pois ~ x,
                                family = "poisson", data = exdata, xcoord = xcoord, ycoord = ycoord,
                                spcov_initial = spcov_initial_known,
                                dispersion_initial = dispersion_initial_known,
                                random = ~group,
                                randcov_initial = randcov_initial_known,
                                anisotropy = TRUE
    ),
    NA
  )
  expect_equal(as.numeric(coef(spmod_with_randcov, type = "spcov")), as.numeric(spcov_initial_known$initial))
  expect_equal(unname(coef(spmod_with_randcov, type = "randcov")[["1 | group"]]), 0.2)
})

test_that("predict() errors informatively when newdata is missing coordinate columns", {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))

  spmod1 <- spglm(abs(y) ~ x, "Gamma", exdata, spcov_type = "exponential", xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml")

  # sanity: predict() still works normally when both coordinate columns are present
  expect_vector(predict(spmod1, newdata = newexdata))

  newexdata_nocoord <- newexdata
  newexdata_nocoord$xcoord <- NULL
  newexdata_nocoord$ycoord <- NULL
  expect_error(predict(spmod1, newdata = newexdata_nocoord), "coordinate column")
})

test_that("predict() errors informatively when newdata has NA in a random intercept grouping column", {
  set.seed(11)
  n <- 60
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  grp <- factor(sample(letters[1:5], n, replace = TRUE))
  p <- plogis(0.3 * x1)
  ybin <- rbinom(n, 1, p)
  d <- data.frame(ybin = ybin, x1 = x1, grp = grp, xcoord = xcoord, ycoord = ycoord)

  gmod <- spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~ (1 | grp))

  newdata_good <- data.frame(
    x1 = c(0.5, -0.5), grp = factor(c("a", "b"), levels = levels(grp)),
    xcoord = c(0.3, 0.6), ycoord = c(0.3, 0.6)
  )
  expect_vector(predict(gmod, newdata_good))

  newdata_na <- newdata_good
  newdata_na$grp[1] <- NA
  expect_error(predict(gmod, newdata_na), "Cannot have NA values in predictors.")
  expect_error(predict(gmod, newdata_na, se.fit = TRUE), "Cannot have NA values in predictors.")
})

test_that("spglm() errors informatively when a formula/random/partition_factor variable is not in data", {
  set.seed(12)
  n <- 30
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  grp <- factor(sample(letters[1:3], n, replace = TRUE))
  ybin <- rbinom(n, 1, plogis(0.3 * x1))
  d <- data.frame(ybin = ybin, x1 = x1, grp = grp, xcoord = xcoord, ycoord = ycoord)

  # a same-named object in the calling environment (but not in data) should
  # not be silently picked up via ordinary formula scoping -- it should error
  not_a_col <- rnorm(n)
  not_a_group <- factor(sample(letters[1:3], n, replace = TRUE))

  expect_error(spglm(ybin ~ not_a_col, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "not_a_col.*not found in data")
  expect_error(spglm(not_a_col ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "not_a_col.*not found in data")
  expect_error(spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~not_a_group), "not_a_group.*not found in data")
  expect_error(spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, partition_factor = ~not_a_group), "not_a_group.*not found in data")

  # sanity: a valid call still works
  expect_s3_class(spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord), "spglm")
})

test_that("spglm() formula supports . as shorthand for all predictors, excluding coordinates", {
  # see the matching test in test-splm.R for the full rationale
  set.seed(13)
  n <- 30
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  ybin <- rbinom(n, 1, plogis(0.3 * x1))
  d <- data.frame(ybin = ybin, x1 = x1, xcoord = xcoord, ycoord = ycoord)

  mod_dot <- spglm(ybin ~ ., family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  mod_explicit <- spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
  expect_equal(names(coef(mod_dot)), names(coef(mod_explicit)))
  expect_equal(unname(coef(mod_dot)), unname(coef(mod_explicit)))
  expect_false(any(c("xcoord", "ycoord") %in% colnames(model.matrix(mod_dot))))

  d_sf <- sf::st_as_sf(d, coords = c("xcoord", "ycoord"))
  mod_sf_dot <- spglm(ybin ~ ., family = "binomial", data = d_sf, spcov_type = "exponential")
  expect_equal(unname(coef(mod_sf_dot)), unname(coef(mod_explicit)))

  expect_error(
    spglm(ybin ~ x1, family = "binomial", data = d, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, random = ~.),
    "not supported in random"
  )
})

test_that("a user-supplied local$index must match the length of the non-missing response", {
  # see the matching test in test-splm.R for the full rationale
  set.seed(12)
  n <- 40
  x1 <- rnorm(n)
  xcoord <- runif(n)
  ycoord <- runif(n)
  ybin <- rbinom(n, 1, plogis(0.3 * x1))
  d <- data.frame(ybin = ybin, x1 = x1, xcoord = xcoord, ycoord = ycoord)

  na_rows <- c(3, 10, 17, 25, 33)
  d_na <- d
  d_na$ybin[na_rows] <- NA
  n_obs <- n - length(na_rows)

  full_index <- sample(1:4, n, replace = TRUE) # wrong length: matches original n, not n_obs
  observed_index <- which(!is.na(d_na$ybin))
  correct_index <- full_index[observed_index] # right length: matches n_obs

  expect_error(
    get_data_object_spglm(
      formula = ybin ~ x1, family = "binomial", data = d_na, spcov_initial = spcov_initial("exponential"),
      xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml", anisotropy = FALSE,
      random = NULL, randcov_initial = NULL, partition_factor = NULL,
      local = list(index = full_index), range_constrain = FALSE
    ),
    "local\\$index must have the same length"
  )
  expect_error(
    spglm(ybin ~ x1, family = "binomial", data = d_na, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = list(index = full_index)),
    "local\\$index must have the same length"
  )

  data_object <- get_data_object_spglm(
    formula = ybin ~ x1, family = "binomial", data = d_na, spcov_initial = spcov_initial("exponential"),
    xcoord = "xcoord", ycoord = "ycoord", estmethod = "reml", anisotropy = FALSE,
    random = NULL, randcov_initial = NULL, partition_factor = NULL,
    local = list(index = correct_index), range_constrain = FALSE
  )
  expect_length(data_object$local_index, n_obs)
  expect_equal(as.vector(data_object$local_index), as.vector(correct_index))

  expect_s3_class(
    spglm(ybin ~ x1, family = "binomial", data = d_na, spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, local = list(index = correct_index)),
    "spglm"
  )
})
