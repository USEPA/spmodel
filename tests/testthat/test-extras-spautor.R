test_that("blank test", {
  expect_null(NULL)
})

test_local <- FALSE # FALSE for CRAN

#### local check
if (test_local) {
  load(file = system.file("extdata", "exdata.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_poly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Mpoly.rda", package = "spmodel"))
  load(file = system.file("extdata", "exdata_Upoly.rda", package = "spmodel"))

  ################################
  ########### row standardized
  ################################

  test_that("the model runs for connected sites", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type), NA)
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, estmethod = "ml"), NA)
    spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
    expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val), NA)
    expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)

    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type), NA)
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, estmethod = "ml"), NA)
    spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
    expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val), NA)
    expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
  })

  test_that("the model runs for connected sites given W", {
    spcov_type <- "car"
    W <- 1 * Matrix::Matrix(sf::st_intersects(exdata_poly, sparse = FALSE), sparse = TRUE)
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = W), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = W, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = W), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = W, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = W), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = W, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = W), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = W, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for connected sites (partition group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, partition_factor = ~group), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, partition_factor = ~group, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, partition_factor = ~group), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, partition_factor = ~group, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for connected sites (random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, estmethod = "ml"), NA)
    }
    spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
    randcov_initial_val <- randcov_initial(group = 1, known = "group")
    expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val), NA)

    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val, estmethod = "ml"), NA)
    }

    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, estmethod = "ml"), NA)
    }

    randcov_initial_val <- randcov_initial(group = 1, known = "group")
    expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val), NA)

    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for connected sites (partition and random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for unconnected sites", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for unconnected sites (partition group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, partition_factor = ~group), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, partition_factor = ~group, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, partition_factor = ~group), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, partition_factor = ~group, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for unconnected sites (random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for unconnected sites (partition and random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for missing sites", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for missing sites (partition group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, partition_factor = ~group, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, partition_factor = ~group, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for missing sites (random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, random = ~group, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, random = ~group, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for connected sites (partition and random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, partition_factor = ~group), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.5, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, estmethod = "ml"), NA)
    }
  })

  ################################
  ########### row unstandardized
  ################################

  test_that("the model runs for connected sites", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for connected sites given W", {
    spcov_type <- "car"
    W <- 1 * Matrix::Matrix(sf::st_intersects(exdata_poly, sparse = FALSE), sparse = TRUE)
    Wp <- W / rowSums(W)
    M <- 1 / rowSums(W)
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = W, row_st = FALSE), NA)
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = Wp, row_st = FALSE, M = M), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = W, row_st = FALSE, estmethod = "ml"), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = Wp, row_st = FALSE, M = M, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = W, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = Wp, row_st = FALSE, M = M), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = W, row_st = FALSE, estmethod = "ml"), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = Wp, row_st = FALSE, M = M, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = W, row_st = FALSE), NA)
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = Wp, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = W, row_st = FALSE, estmethod = "ml"), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, W = Wp, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = W, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = Wp, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = W, row_st = FALSE, estmethod = "ml"), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, W = Wp, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for connected sites (partition group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for connected sites (random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
    spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
    randcov_initial_val <- randcov_initial(group = 1, known = "group")
    expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val, row_st = FALSE), NA)

    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
    spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
    randcov_initial_val <- randcov_initial(group = 1, known = "group")
    expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, randcov_initial = randcov_initial_val, row_st = FALSE), NA)
  })

  test_that("the model runs for connected sites (partition and random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }

    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_poly, spcov_type, random = ~group, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for unconnected sites", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, row_st = FALSE, estmethod = "reml"), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, row_st = FALSE, estmethod = "reml"), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for unconnected sites (partition group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE, estmethod = "reml"), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE, estmethod = "reml"), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for unconnected sites (random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE, estmethod = "reml"), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE, estmethod = "reml"), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for unconnected sites (partition and random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_type, random = ~group, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_Upoly, spcov_initial = spcov_initial_val, random = ~group, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for missing sites", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for missing sites (partition group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, partition_factor = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, partition_factor = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("the model runs for missing sites (random group)", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }


    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, row_st = FALSE), NA)
    if (test_local) {
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_type, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
      spcov_initial_val <- spcov_initial(spcov_type, de = 1, ie = 1, range = 0.1, known = "de")
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE), NA)
      expect_error(spautor(y ~ x, exdata_Mpoly, spcov_initial = spcov_initial_val, random = ~group, row_st = FALSE, estmethod = "ml"), NA)
    }
  })

  test_that("examples run", {
    spmod <- spautor(y ~ x, exdata_poly, "car")
    expect_output(print(spmod))
    expect_error(summary(spmod), NA)
    expect_error(tidy(spmod), NA)
    # random effects
    spmod <- expect_error(spautor(y ~ x, exdata_poly, "car", random = ~group), NA)
    spmod <- expect_error(spautor(y ~ x, exdata_poly, "car", random = ~ (x | group) + (x | subgroup)), NA)
    # partition factor
    spmod <- expect_error(spautor(y ~ x, exdata_poly, "car", partition_factor = ~group), NA)
    # combining
    spmod <- expect_error(spautor(y ~ x, exdata_poly, "car", random = ~group, partition_factor = ~group), NA)
    # spcov_initial
    spcov_initial_val <- spcov_initial("car", range = 0.5, known = "range")
    spmod <- expect_error(spautor(y ~ x, exdata_poly, spcov_initial = spcov_initial_val), NA)
    # randcov_initial ("group" is shorthand for "1 | group")
    randcov_initial_val <- randcov_initial(1, nm = "group", known = "group")
    spmod <- expect_error(spautor(y ~ x, exdata_poly, "car", randcov_initial = randcov_initial_val), NA)
  })

  test_that("errors occur", {
    W <- 1 * sf::st_intersects(exdata_poly, sparse = FALSE)
    diag(W) <- 0
    M <- seq_len(NROW(exdata_poly))
    expect_error(spautor(y ~ x, exdata_poly, "car", W = W, row_st = FALSE, M = M))
    expect_error(spautor(y ~ x, exdata_poly, "car", partition_factor = ~ group + subgroup))
    expect_error(suppressWarnings(spautor(y ~ as.factor(x) + group, exdata_poly, "car")))
    exdata_poly3 <- exdata_poly
    exdata_poly3$y <- as.character(exdata_poly3$y)
    expect_error(spautor(y ~ x, exdata_poly3, "car"))
    expect_error(spautor(as.character(y) ~ x, exdata_poly, "car"))
    exdata_poly3$y <- as.factor(exdata_poly3$y)
    expect_error(spautor(y ~ x, exdata_poly3, "car"))
    expect_error(spautor(as.factor(y) ~ x, exdata_poly, "car"))
    expect_error(spautor(y ~ x, exdata_poly, "xyz"))
    expect_error(spautor(y ~ x, exdata_poly, "exponential"))
    expect_error(spautor(y ~ x, exdata, "car"))
    expect_error(spautor(y ~ x, exdata_poly, "car", estmethod = "xyz"))
    exdata_poly4 <- exdata_poly
    exdata_poly4$x2 <- exdata_poly4$x
    expect_error(suppressWarnings(spautor(y ~ x + x2, exdata_poly4, "car")))
    exdata_poly4$x[1] <- NA
    expect_error(spautor(y ~ x, exdata_poly4, "car"))
  })

  test_that("additional covr", {
    # spatial profiling
    spcov_initial_val <- spcov_initial("car", ie = 0.5) # just can't be zero and known
    expect_error(spautor(formula = y ~ x, data = exdata_poly, spcov_initial = spcov_initial_val), NA)

    # partition factor with one level and random effects
    exdata_poly$part <- as.factor(1)
    expect_error(spautor(formula = y ~ x, data = exdata_poly, "car", random = ~group, partition_factor = ~part), NA)
    exdata_Upoly$part <- as.factor(1)
    expect_error(spautor(formula = y ~ x, data = exdata_Upoly, "car", random = ~group, partition_factor = ~part), NA)
    exdata_Mpoly$part <- as.factor(1)
    expect_error(spautor(formula = y ~ x, data = exdata_Mpoly, "car", random = ~group, partition_factor = ~part), NA)
  })

  test_that("messages occur", {
    expect_message(spautor(y ~ x, exdata_poly))
    spcov_initial_val <- spcov_initial("car")
    expect_message(spautor(y ~ x, exdata_poly, "car", spcov_initial = spcov_initial_val))
  })

  test_that("no variance error works", {
    exdata_poly$novar <- 1
    expect_error(spautor(novar ~ x, exdata_poly, "car"))
  })

  test_that("range negative works", {
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, range_positive = FALSE), NA)
    spcov_type <- "sar"
    expect_error(spautor(y ~ x, exdata_poly, spcov_type, range_positive = FALSE), NA)
  })

  test_that("emmeans works", {
    spcov_type <- "car"
    spmod <- spautor(y ~ x * group, exdata_poly, spcov_type = spcov_type, estmethod = "reml")
    expect_equal(as.matrix(model.frame(delete.response(terms(spmod)), spmod$data[spmod$observed_index, , drop = FALSE])), as.matrix(emmeans::recover_data(spmod)))
    expect_error(emmeans::emmeans(spmod, ~ group, by = "x"), NA)
  })

  test_that("emmeans works missing", {
    spcov_type <- "car"
    spmod <- spautor(y ~ x * group, exdata_Mpoly, spcov_type = spcov_type, estmethod = "reml")
    expect_equal(as.matrix(model.frame(delete.response(terms(spmod)), spmod$data[spmod$observed_index, , drop = FALSE])), as.matrix(emmeans::recover_data(spmod)))
    expect_error(emmeans::emmeans(spmod, ~ group, by = "x"), NA)
  })

  test_that("point distance works missing", {
    exdata_sf <- st_as_sf(exdata, coords = c("xcoord", "ycoord"), crs = NA)
    spcov_type <- "car"
    expect_error(spautor(y ~ x, exdata_sf, spcov_type = spcov_type, cutoff = 1), NA)
    expect_error(spautor(y ~ x, exdata_sf, spcov_type = spcov_type, cutoff = 1, row_st = FALSE), NA)
    expect_error(spautor(y ~ x, exdata_sf, spcov_type = spcov_type, cutoff = NULL)) # can't be NULL
    expect_error(spautor(y ~ x, exdata_sf, spcov_type = spcov_type, cutoff = 1e-8)) # too small of distance so no neighbors
  })
}
