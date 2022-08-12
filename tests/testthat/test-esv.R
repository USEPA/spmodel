load(file = system.file("extdata", "exdata.rda", package = "spmodel"))

test_that("esv works", {

  # regular implementation
  esv1 <- esv(y ~ x, exdata, xcoord, ycoord)
  expect_s3_class(esv1, "data.frame")
  expect_equal(NROW(esv1), 15)
  expect_equal(NCOL(esv1), 4)

  esv1_q <- esv(y ~ x, exdata, "xcoord", "ycoord")
  expect_s3_class(esv1_q, "data.frame")
  expect_equal(NROW(esv1_q), 15)
  expect_equal(NCOL(esv1_q), 4)

  # specifying bins and cutoff
  esv2 <- esv(y ~ x, exdata, xcoord, ycoord, bins = 30, cutoff = 5)
  expect_s3_class(esv2, "data.frame")
  expect_equal(NROW(esv2), 30)
  expect_equal(NCOL(esv1), 4)
  dist_matrix <- spdist(exdata, "xcoord", "ycoord")

  # specifying distance matrix
  esv3 <- esv(y ~ x, exdata, dist_matrix = dist_matrix)
  expect_s3_class(esv3, "data.frame")
  expect_equal(NROW(esv3), 15)
  expect_equal(NCOL(esv1), 4)

  # specifying partition factor
  esv4 <- esv(y ~ x, exdata, xcoord, ycoord, partition_factor = ~group)
  expect_s3_class(esv1, "data.frame")
  expect_equal(NROW(esv1), 15)
  expect_equal(NCOL(esv1), 4)
  expect_false(identical(esv1, esv4)) # make sure results are not identical to full esv

  # quoting works
  esv1_2 <- esv(y ~ x, exdata, "xcoord", "ycoord")
  expect_equal(esv1, esv1_2)

  # works with sf object
  exdata_sf <- sf::st_as_sf(exdata, coords = c("xcoord", "ycoord"))
  expect_error(esv(y ~ x, exdata_sf), NA)
})
