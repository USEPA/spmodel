test_that("blank test", {
  expect_null(NULL)
})

test_local <- FALSE # FALSE for CRAN

if (test_local) {
  set.seed(1)

  # SPMODEL PACKAGE NEEDS TO BE INSTALLED VIA DEVTOOLS::INSTALL() BEFORE RUNNING TESTS IF THOSE TESTS HAVE PARALLELIZATION

  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

  # add variables
  n <- NROW(exdata_poly)
  exdata_poly$bern <- rbinom(n, size = 1, prob = 0.5)
  exdata_poly$bernfac <- factor(ifelse(exdata_poly$bern == 0, "a", "b"))
  exdata_poly$size <- 5
  exdata_poly$bin <- rbinom(n, size = exdata_poly$size, prob = 0.5)
  exdata_poly$prop <- runif(n)
  exdata_poly$count <- rpois(n, lambda = 5)
  exdata_poly$cont <- rgamma(n, shape = 1, rate = 1)
  exdata_poly$offset <- 2

  # save W and M
  W <- sf::st_intersects(exdata_poly, sparse = FALSE)
  diag(W) <- 0
  W <- 1 * Matrix::Matrix(W, sparse = TRUE)
  W_rowsums <- Matrix::rowSums(W)
  Wbar <- W / W_rowsums
  M <- Matrix(diag(NROW(W))) # for row-unstandardized / as matrix
  Mvec <- 1 / W_rowsums # for row-standardized / as vector

  test_that("the model runs for binomial data", {
    spgmod <- spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml")
    expect_s3_class(spgmod, "spgautor")
    expect_vector(AUROC(spgmod))
    expect_error(spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml", range_positive = FALSE), NA)
    expect_error(spgautor(bernfac ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "ml"), NA)
    expect_error(spgautor(cbind(bin, size) ~ x, family = "binomial", data = exdata_poly, spcov_type = "sar", estmethod = "reml"), NA)
    # causes an error with ml estimation as de is near zero and ie is zero, which makes inverse unstable
    # need to implement a diagonal tolerance for gautor models
    # expect_error(spgautor(cbind(bin, size) ~ x, family = "binomial", data = exdata_poly, spcov_type = "sar", estmethod = "reml"), NA)
    spgmod <- spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml")
    expect_vector(AUROC(spgmod))

    # complicated models
    expect_error(spgautor(bern ~ x,
      family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml",
      W = W, row_st = FALSE, M = M
    ), NA)
    expect_error(spgautor(bernfac ~ x,
      family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml",
      W = W, row_st = FALSE, M = M
    ), NA)
    expect_error(spgautor(cbind(bin, size) ~ x + offset(offset),
      family = "binomial", data = exdata_poly, spcov_type = "sar", estmethod = "ml",
      random = ~group
    ), NA)
    expect_error(spgautor(y > 0 ~ x,
      family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml",
      W = W, row_st = FALSE, M = M
    ), NA)
  })

  test_that("the model runs for proportion data", {
    expect_error(spgautor(prop ~ x, family = "beta", data = exdata_poly, spcov_type = "car", estmethod = "reml"), NA)
    expect_error(spgautor(prop ~ x, family = "beta", data = exdata_poly, spcov_type = "car", estmethod = "reml", range_positive = FALSE), NA)
    expect_error(spgautor(prop ~ x, family = beta, data = exdata_poly, spcov_type = "sar", estmethod = "reml"), NA)

    # complicated models
    expect_error(spgautor(prop ~ x,
      family = "beta", data = exdata_poly, spcov_type = "car", estmethod = "reml",
      random = ~ group + subgroup
    ), NA)
    expect_error(spgautor(prop ~ x + offset(offset),
      family = beta, data = exdata_poly, spcov_type = "sar", estmethod = "reml",
      W = W, partition_factor = ~group
    ), NA)
  })

  test_that("the model runs for count data", {
    spgmod <- spgautor(count ~ x, family = poisson, data = exdata_poly, spcov_type = "sar", estmethod = "reml")
    expect_s3_class(spgmod, "spgautor")
    expect_error(AUROC(spgmod))
    expect_error(spgautor(count ~ x, family = "nbinomial", data = exdata_poly, spcov_type = "car", estmethod = "ml"), NA)
    expect_error(spgautor(count ~ x, family = "nbinomial", data = exdata_poly, spcov_type = "car", estmethod = "ml", range_positive = FALSE), NA)

    # complicated models
    expect_error(spgautor(count ~ x,
      family = poisson, data = exdata_poly, spcov_type = "sar", estmethod = "reml",
      random = ~ (x | subgroup)
    ), NA)
    expect_error(spgautor(count ~ x + offset(offset),
      family = "nbinomial", data = exdata_poly, spcov_type = "car", estmethod = "ml",
      W = Wbar, row_st = FALSE, M = Mvec
    ), NA)
  })

  test_that("the model runs for continuous data", {
    expect_error(spgautor(cont ~ x, family = "Gamma", data = exdata_poly, spcov_type = "car", estmethod = "reml"), NA)
    expect_error(spgautor(cont ~ x, family = inverse.gaussian, data = exdata_poly, spcov_type = "sar", estmethod = "ml"), NA)
    expect_error(spgautor(cont ~ x, family = inverse.gaussian, data = exdata_poly, spcov_type = "car", estmethod = "reml", range_positive = FALSE), NA)
    # SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
    expect_error(spgautor(cont ~ x, family = gaussian, data = exdata_poly, spcov_type = "car", estmethod = "reml"))

    # complicated models
    expect_error(spgautor(cont ~ x + offset(offset),
      family = "Gamma", data = exdata_poly, spcov_type = "car", estmethod = "reml",
      W = Wbar, row_st = FALSE, M = Mvec
    ), NA)
    expect_error(spgautor(cont ~ x,
      family = inverse.gaussian, data = exdata_poly, spcov_type = "sar", estmethod = "ml",
      random = ~subgroup, partition_factor = ~group
    ), NA)
    ## SHOULD BE AN ERROR AS GAUSSIAN FAMILY REMOVED
    expect_error(spgautor(cont ~ x,
      family = gaussian, data = exdata_poly, spcov_type = "car", estmethod = "reml",
      random = ~subgroup
    ))
  })

  test_that("the model runs on other data sets", {
    expect_error(spgautor(abs(y) ~ x, family = "Gamma", data = exdata_Mpoly, spcov_type = "car", estmethod = "reml"), NA)
    expect_error(spgautor(abs(y) ~ x, family = Gamma, data = exdata_Upoly, spcov_type = "sar", estmethod = "reml"), NA)

    # complicated models
    expect_error(spgautor(abs(y) ~ x,
      family = "Gamma", data = exdata_Mpoly, spcov_type = "car", estmethod = "reml",
      W = W, partition_factor = ~subgroup
    ), NA)

    W <- sf::st_intersects(exdata_Upoly, sparse = FALSE)
    diag(W) <- 0
    W <- 1 * Matrix::Matrix(W, sparse = TRUE)
    expect_error(spgautor(abs(y) ~ x,
      family = Gamma, data = exdata_Upoly, spcov_type = "sar", estmethod = "reml",
      W = W, random = ~subgroup
    ), NA)
  })

  test_that("emmeans works", {
    spcov_type <- "car"
    spgmod <- spgautor(abs(y) ~ x * group, family = "Gamma", exdata_poly, spcov_type = spcov_type, estmethod = "reml")
    expect_equal(as.matrix(model.frame(delete.response(terms(spgmod)), spgmod$data[spgmod$observed_index, , drop = FALSE])), as.matrix(emmeans::recover_data(spgmod)))
    expect_error(emmeans::emmeans(spgmod, ~ group, by = "x"), NA)
  })

  test_that("emmeans works missing", {
    spcov_type <- "car"
    spgmod <- spgautor(abs(y) ~ x * group, family = "Gamma", exdata_Mpoly, spcov_type = spcov_type, estmethod = "reml")
    expect_equal(as.matrix(model.frame(delete.response(terms(spgmod)), spgmod$data[spgmod$observed_index, , drop = FALSE])), as.matrix(emmeans::recover_data(spgmod)))
    expect_error(emmeans::emmeans(spgmod, ~ group, by = "x"), NA)
  })

  test_that("point distance works missing", {
    exdata_sf <- st_as_sf(exdata, coords = c("xcoord", "ycoord"), crs = NA)
    spcov_type <- "sar"
    expect_error(spgautor(abs(y) ~ x, family = "Gamma", exdata_sf, spcov_type = spcov_type, cutoff = 1), NA)
    expect_error(spgautor(abs(y) ~ x, family = Gamma, exdata_sf, spcov_type = spcov_type, cutoff = 1, row_st = FALSE), NA)
    expect_error(spgautor(abs(y) ~ x, family = "Gamma", exdata_sf, spcov_type = spcov_type, cutoff = NULL)) # can't be NULL
    expect_error(spgautor(abs(y) ~ x, family = Gamma, exdata_sf, spcov_type = spcov_type, cutoff = 1e-8)) # too small of distance so no neighbors
  })
}
