load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
load(file = system.file("extdata", "newexdata.rda", package = "spmodel"))
load(system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))

test_that("augment works geo", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  expect_error(augment(spmod), NA)
  expect_error(augment(spmod, se_fit = TRUE), NA)
  expect_s3_class(augment(spmod), "tbl")
  expect_error(augment(spmod, newdata = newexdata), NA)
  augment(spmod, newdata = newexdata, interval = "confidence")
  augment(spmod, newdata = newexdata, interval = "prediction")
  expect_error(augment(spmod, newdata = newexdata, se_fit = TRUE), NA)
  augment(spmod, newdata = newexdata, interval = "confidence", se_fit = TRUE)
  augment(spmod, newdata = newexdata, interval = "prediction", se_fit = TRUE)
  expect_s3_class(augment(spmod, newdata = newexdata), "tbl")
})

test_that("augment works auto", {
  spmod <- spautor(y ~ x, exdata_Mpoly, "car")
  expect_error(augment(spmod), NA)
  expect_error(augment(spmod, se_fit = TRUE), NA)
  expect_s3_class(augment(spmod), "tbl")
  expect_s3_class(augment(spmod), "sf")
  expect_error(augment(spmod, newdata = spmod$newdata), NA)
  augment(spmod, newdata = spmod$newdata, interval = "confidence")
  augment(spmod, newdata = spmod$newdata, interval = "prediction")
  expect_error(augment(spmod, newdata = spmod$newdata, se_fit = TRUE), NA)
  augment(spmod, newdata = spmod$newdata, interval = "confidence", se_fit = TRUE)
  augment(spmod, newdata = spmod$newdata, interval = "prediction", se_fit = TRUE)
  expect_s3_class(augment(spmod, newdata = spmod$newdata), "tbl")
  expect_s3_class(augment(spmod, newdata = spmod$newdata), "sf")
})

test_that("augment works with drop", {
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  aug_mod <- augment(spmod) # default drop = TRUE
  expect_equal(NCOL(aug_mod), 9)
  aug_mod <- augment(spmod, drop = FALSE)
  expect_equal(NCOL(aug_mod), 11)

  aug_pred <- augment(spmod, newdata = newexdata)
  expect_equal(NCOL(aug_pred), 6) # default drop = FALSE

  spmod <- spautor(y ~ x, exdata_Mpoly, "car")
  aug_mod <- augment(spmod) # default drop = TRUE
  expect_equal(NCOL(aug_mod), 8)
  aug_mod <- augment(spmod, drop = FALSE)
  expect_equal(NCOL(aug_mod), 11)

  aug_pred <- augment(spmod, newdata = spmod$newdata)
  expect_equal(NCOL(aug_pred), 7) # default drop = FALSE
})

test_that("augment works with types", {
  exdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"))
  newexdata_sf <- sf::st_as_sf(newexdata, coords = c("xcoord", "ycoord"))

  # newdata output type same as data output type

  # df fit df pred
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  aug_mod <- augment(spmod) # default drop = TRUE
  expect_true(inherits(aug_mod, "tbl"))
  expect_false(inherits(aug_mod, "sf"))
  aug_mod <- augment(spmod, drop = FALSE)
  expect_true(inherits(aug_mod, "tbl"))
  expect_false(inherits(aug_mod, "sf"))
  aug_pred <- augment(spmod, newdata = newexdata)
  expect_true(inherits(aug_pred, "tbl"))
  expect_false(inherits(aug_pred, "sf"))

  # sf fit sf pred
  spmod <- splm(y ~ x, exdata_sf, "exponential", xcoord, ycoord)
  aug_mod <- augment(spmod) # default drop = TRUE
  expect_true(inherits(aug_mod, "tbl"))
  expect_true(inherits(aug_mod, "sf"))
  aug_mod <- augment(spmod, drop = FALSE)
  expect_true(inherits(aug_mod, "tbl"))
  expect_true(inherits(aug_mod, "sf"))
  aug_pred <- augment(spmod, newdata = newexdata_sf)
  expect_true(inherits(aug_pred, "tbl"))
  expect_true(inherits(aug_pred, "sf"))

  # df fit sf pred
  spmod <- splm(y ~ x, exdata, "exponential", xcoord, ycoord)
  aug_mod <- augment(spmod) # default drop = TRUE
  expect_true(inherits(aug_mod, "tbl"))
  expect_false(inherits(aug_mod, "sf"))
  aug_mod <- augment(spmod, drop = FALSE)
  expect_true(inherits(aug_mod, "tbl"))
  expect_false(inherits(aug_mod, "sf"))
  aug_pred <- augment(spmod, newdata = newexdata_sf)
  expect_true(inherits(aug_pred, "tbl"))
  expect_false(inherits(aug_pred, "sf"))

  # sf fit df pred
  newexdata$.xcoord <- newexdata$xcoord
  newexdata$.ycoord <- newexdata$ycoord
  spmod <- splm(y ~ x, exdata_sf, "exponential", xcoord, ycoord)
  aug_mod <- augment(spmod) # default drop = TRUE
  expect_true(inherits(aug_mod, "tbl"))
  expect_true(inherits(aug_mod, "sf"))
  aug_mod <- augment(spmod, drop = FALSE)
  expect_true(inherits(aug_mod, "tbl"))
  expect_true(inherits(aug_mod, "sf"))
  aug_pred <- augment(spmod, newdata = newexdata)
  expect_true(inherits(aug_pred, "tbl"))
  expect_true(inherits(aug_pred, "sf"))



  # sf in
  spmod <- spautor(y ~ x, exdata_Mpoly, "car")
  expect_true(inherits(aug_mod, "tbl"))
  expect_true(inherits(aug_mod, "sf"))
  aug_mod <- augment(spmod, drop = FALSE)
  expect_true(inherits(aug_mod, "tbl"))
  expect_true(inherits(aug_mod, "sf"))
  aug_pred <- augment(spmod, newdata = spmod$newexdata)
  expect_true(inherits(aug_pred, "tbl"))
  expect_true(inherits(aug_pred, "sf"))

  # df in
  W <- 1 * sf::st_intersects(exdata_Mpoly, sparse = FALSE)
  diag(W) <- 0
  exdata_Mpoly_df <- st_drop_geometry(exdata_Mpoly)
  spmod <- spautor(y ~ x, exdata_Mpoly_df, "car", W = W)
  aug_mod <- augment(spmod)
  expect_true(inherits(aug_mod, "tbl"))
  expect_false(inherits(aug_mod, "sf"))
  aug_mod <- augment(spmod, drop = FALSE)
  expect_true(inherits(aug_mod, "tbl"))
  expect_false(inherits(aug_mod, "sf"))
  aug_pred <- augment(spmod, newdata = spmod$newdata)
  expect_true(inherits(aug_pred, "tbl"))
  expect_false(inherits(aug_pred, "sf"))
})
