# Snapshot tests pinning the optimizer *starting values* that
# cov_initial_search()/cov_initial_search_glm() choose, one representative
# case per covariance-type group/estmethod/random-effect combination.

capture_cov_initial <- function(fn_name, expr) {
  # trace()'s exit hook is evaluated in the traced function's own call frame,
  # whose lexical scope is the spmodel namespace -- not this helper's local
  # frame or test_that()'s
  suppressMessages(trace(fn_name,
    exit = quote(assign(".captured_cov_initial", returnValue(), envir = globalenv())),
    print = FALSE, where = asNamespace("spmodel")
  ))
  on.exit(suppressMessages(untrace(fn_name, where = asNamespace("spmodel"))))
  force(expr)
  captured <- get(".captured_cov_initial", envir = globalenv())
  rm(".captured_cov_initial", envir = globalenv())
  captured
}

load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))

set.seed(1)
exdata$bern <- rbinom(NROW(exdata), size = 1, prob = 0.5)
exdata_poly$bern <- rbinom(NROW(exdata_poly), size = 1, prob = 0.5)

test_that("splm starting grid is unchanged (exponential)", {
  captured <- capture_cov_initial(
    "cov_initial_search",
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search",
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "sv-wls")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search",
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "sv-cl")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search",
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$randcov_initial_val$initial))
})

test_that("splm starting grid is unchanged (anisotropy)", {
  # exercises the rotate/pi-rotate ambiguity resolution inside eval_grid(),
  # not covered by any of the (non-anisotropic) blocks above
  captured <- capture_cov_initial(
    "cov_initial_search",
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", anisotropy = TRUE)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
})

test_that("splm starting grid is unchanged (none/ie)", {
  captured <- capture_cov_initial(
    "cov_initial_search",
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search",
    suppressWarnings(splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml", random = ~group)) # vcov_theta pd warning
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$randcov_initial_val$initial))
})

test_that("splm starting grid is unchanged (matern)", {
  captured <- capture_cov_initial(
    "cov_initial_search",
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search",
    splm(y ~ x, exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml", random = ~group)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$randcov_initial_val$initial))
})

test_that("spautor starting grid is unchanged (car/sar)", {
  captured <- capture_cov_initial(
    "cov_initial_search",
    spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search",
    spautor(y ~ x, exdata_poly, spcov_type = "car", estmethod = "reml", random = ~group)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$randcov_initial_val$initial))
})

test_that("spglm starting grid is unchanged (exponential)", {
  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$dispersion_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", random = ~group)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$randcov_initial_val$initial))
})

test_that("spglm starting grid is unchanged (anisotropy)", {
  # exercises the rotate/pi-rotate ambiguity resolution inside
  # eval_grid_glm(), not covered by any of the (non-anisotropic) blocks above
  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "exponential", estmethod = "reml", anisotropy = TRUE)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
})

test_that("spglm starting grid is unchanged (none/ie)", {
  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "none", estmethod = "reml", random = ~group)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$randcov_initial_val$initial))
})

test_that("spglm starting grid is unchanged (matern)", {
  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spglm(bern ~ x, family = binomial, data = exdata, xcoord = xcoord, ycoord = ycoord, spcov_type = "matern", estmethod = "reml", random = ~group)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$randcov_initial_val$initial))
})

test_that("spgautor starting grid is unchanged (car/sar)", {
  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml")
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))

  captured <- capture_cov_initial(
    "cov_initial_search_glm",
    spgautor(bern ~ x, family = binomial, data = exdata_poly, spcov_type = "car", estmethod = "reml", random = ~group)
  )
  expect_snapshot(print(captured$spcov_initial_val$initial))
  expect_snapshot(print(captured$randcov_initial_val$initial))
})
