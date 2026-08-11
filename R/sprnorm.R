#' Simulate a spatial normal (Gaussian) random variable
#'
#' @description Simulate a spatial normal (Gaussian) random variable with a specific
#'   mean and covariance structure.
#'
#' @param spcov_params An [spcov_params()] object.
#' @param mean A numeric vector representing the mean. \code{mean} must have length 1
#'   (in which case it is recycled) or length equal
#'   to the number of rows in \code{data}. The default is \code{0}.
#' @param samples The number of independent samples to generate. The default
#'   is \code{1}.
#' @param data A data frame or \code{sf} object containing spatial information.
#' @param randcov_params A [randcov_params()] object.
#' @param partition_factor A formula indicating the partition factor.
#' @param ... Other arguments. Not used (needed for generic consistency).
#' @param xcoord Name of the column in \code{data} representing the x-coordinate.
#'   Can be quoted or unquoted. Not required if \code{data} are an \code{sf}
#'   object.
#' @param ycoord Name of the column in \code{data} representing the y-coordinate.
#'   Can be quoted or unquoted. Not required if \code{data} are an \code{sf}
#'   object.
#' @param local An optional logical or list controlling the big data approximation.
#'   If omitted, \code{local} is set
#'   to \code{TRUE} or \code{FALSE} based on the desired sample size (the number of
#'   non-missing observations in \code{data}) -- if the desired sample size exceeds 5,000,
#'   \code{local} is set to \code{TRUE}. Otherwise it is set to \code{FALSE}.
#'   \code{local} is also set to \code{FALSE} when \code{spcov_type} is \code{"none"}
#'   and there are no random effects specified via \code{random}.
#'   If \code{FALSE}, no big data approximation is implemented.
#'   If a list is provided, \code{local$approximation} selects which big data
#'   approximation is used and can take on the values
#'   \itemize{
#'     \item \code{"low-rank"}: a base sample is drawn from the data, the
#'       remaining locations are split into blocks, and each block is
#'       simulated conditional on the base sample alone. Blocks are
#'       assumed conditionally independent of one another given the base
#'       sample.
#'       \itemize{
#'         \item \code{method_base}: Whether the data simulated is restricted
#'           to a base sample at all. If \code{method_base = "all"}, no big
#'           data approximation is applied (equivalent to \code{local = FALSE});
#'           this happens automatically whenever \code{size_base} is at least
#'           the desired sample size. If \code{method_base = "base"}, the data
#'           is subset to \code{size_base} locations (ordered via
#'           \code{reorder_base}) to form the base sample. The default is
#'           \code{"base"}.
#'         \item \code{reorder_base}: The data reordering approach used prior to splitting
#'           into base and new sets. If \code{reorder_base = "none"}, no reordering
#'           is applied to the data. If \code{reorder_base = "random"}, the data order is
#'           randomly reshuffled. If \code{reorder_base = "grts"}, the data order is
#'           randomly generated using the GRTS algorithm for spatially balanced
#'           sampling via \code{spsurvey::grts()}. The default is \code{"grts"}.
#'         \item \code{size_base}: The number of data observations used for the base sample.
#'           The default is 5,000. See Details for more.
#'         \item \code{kmeans_new}: For observations outside the base sample, whether
#'           they should be assigned to blocks based on k-means clustering
#'           on the coordinates, with clusters of size approximately equal to
#'           \code{size_new}. The default is \code{FALSE} when \code{reorder_base = "none"}
#'           and \code{TRUE} otherwise.
#'         \item \code{size_new}: The (approximate) number of observations used
#'           for each block. The default is 1,000. See Details for more.
#'         \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'           parallel package is automatically used. The default is \code{FALSE}.
#'         \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'           parallelize over. The default is the number of available cores on your machine.
#'       }
#'       If \code{local$approximation} is \code{"low-rank"} (either explicitly or via
#'       \code{local = TRUE}), defaults for the remaining \code{"low-rank"}
#'       settings are chosen such that \code{local} is transformed into
#'       \code{list(approximation = "low-rank", method_base = "base", size_base = 5000,
#'       reorder_base = "grts", size_new = 1000, kmeans_new = TRUE, parallel = FALSE)}.
#'     \item \code{"vecchia"}: every location is simulated one at a time (in
#'       some order), each conditional on \strong{all} already-simulated
#'       locations (not a single shared base sample. Locations are
#'       never assumed conditionally independent of one another. This is
#'       exact (matches \code{local = FALSE}) when \code{method = "all"};
#'       \code{method = "distance"}/\code{"covariance"} truncate the conditioning set to a fixed number of neighbors
#'       sorted by distance or covariance with the new observation. No parallelization
#'       exists because the algorithm is inherently sequential, as each new observation
#'       depends on previous ones.
#'       \itemize{
#'         \item \code{method}: The neighbor-selection rule used to build each
#'           location's conditioning set once it exceeds \code{size}
#'           already-simulated candidates. Values are \code{"all"} (no truncation,
#'           exact), \code{"distance"} (the \code{size} nearest candidates),
#'           or \code{"covariance"} (the \code{size} candidates with the
#'           highest covariance, in absolute value, with the location being
#'           simulated). Same convention as \code{predict()}'s own
#'           \code{local$method}. The default is \code{"covariance"}. \code{method = "all"} is very computationally
#'           intensive and \code{local = FALSE} should almost always be used instead. 
#'           (\code{method = "all"} primarily exists for numerical verification).
#'         \item \code{size}: The number of neighbors used when \code{method}
#'           is \code{"distance"} or \code{"covariance"}. The default is 30.
#'         \item \code{ordering}: The order locations are simulated in --
#'           \code{"maxmin"}, \code{"middleout"}, \code{"outsidein"},
#'           \code{"coordinate"}, \code{"grts"}, \code{"random"}, or
#'           \code{"none"} (same options as \code{decorrelate()}'s
#'           \code{ordering} argument). The default is \code{"maxmin"}.
#'       }
#'       \code{parallel}/\code{ncores} are not used when \code{local$approximation} is
#'       \code{"vecchia"}.
#'   }
#' 
#'       When \code{local = TRUE}, \code{local} is transformed into
#'       \code{list(approximation = "low-rank", method_base = "base", size_base = 5000,
#'       reorder_base = "grts", method_new = "base", size_new = 1000,
#'       reorder_new = "random", kmeans_new = TRUE, parallel = FALSE)}.
#'       When \code{local} is a list, at least one list element must be provided to
#'       initialize default arguments for the other list elements. See Details for more.
#' @param W Weight matrix specifying the neighboring structure used for car and
#'   sar models. Not required if \code{data} are an \code{sf}
#'   polygon object and \code{W} should be calculated internally (using queen contiguity).
#' @param row_st A logical indicating whether row standardization be performed on
#'   \code{W}. The default is \code{TRUE}.
#' @param M M matrix satisfying the car symmetry condition. The car
#'   symmetry condition states that \eqn{(I - range * W)^{-1}M} is symmetric, where
#'   \eqn{I} is an identity matrix, \eqn{range} is a constant that controls the
#'   spatial dependence, \code{W} is the weights matrix,
#'   and \eqn{^{-1}} represents the inverse operator.
#'   \code{M} is required for car models
#'   when \code{W} is provided and \code{row_st} is \code{FALSE}.  When \code{M},
#'   is required, the default is the identity matrix.
#'
#' @details Random variables are simulated via the product of the covariance matrix's
#'   square (Cholesky) root and independent standard normal random variables
#'   with mean 0 and variance 1, with big data approximations available. It is
#'   nearly the sample computational cost to call \code{sprnorm()} for any value
#'   of \code{samples}.
#'
#'   Only methods for the \code{exponential} and \code{car}
#'   covariance functions are documented here,
#'   but methods exist for all other spatial covariance functions defined in
#'   [spcov_initial()]. Syntax for the \code{exponential} method is the same
#'   as syntax for \code{ie}, \code{spherical}, \code{gaussian}, \code{triangular},
#'   \code{circular}, \code{cubic}, \code{pentaspherical}, \code{cosine}, \code{wave},
#'   \code{jbessel}, \code{gravity}, \code{rquad}, \code{magnetic}, \code{matern},
#'   \code{cauchy}, and \code{pexponential} methods. Syntax for
#'   the \code{car} method is the same as syntax for the \code{sar} method. The
#'   \code{extra} parameter for car and sar models is ignored when all observations have
#'   neighbors.
#'
#'   \code{local} Details: When \code{local$approximation} is \code{"low-rank"}, the big
#'   data approximation works by assigning \code{size_base}
#'   observations to a base sample and then simulating data for the base sample.
#'   The remaining observations are assigned to blocks. For each block, data
#'   are simulated from the conditional distribution given the base sample.
#'   Observations from the same block share conditional covariance while
#'   observations from distinct blocks are assumed conditionally independent
#'   (given the base sample). Parallelization generally further speeds up
#'   computations. When \code{local$approximation} is \code{"vecchia"}, no such
#'   independence assumption is made and each new observation is simulated by
#'   conditioning on prior observations.
#'
#' @return If \code{samples} is 1, a vector of random variables for each row of \code{data}
#'   is returned. If \code{samples} is greater than one, a matrix of random variables
#'   is returned, where the rows correspond to each row of \code{data} and the columns
#'   correspond to independent samples.
#'
#' @export
#'
#' @examples
#' spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
#' sprnorm(spcov_params_val, data = caribou, xcoord = x, ycoord = y)
#' sprnorm(spcov_params_val, mean = 1:30, samples = 5, data = caribou, xcoord = x, ycoord = y)
sprnorm <- function(spcov_params, mean = 0, samples = 1, data, randcov_params, partition_factor,  ...) {
  UseMethod("sprnorm", spcov_params)
}
#' @rdname sprnorm
#' @method sprnorm exponential
#' @export
sprnorm.exponential <- function(spcov_params, mean = 0, samples = 1, data, randcov_params, partition_factor, xcoord, ycoord, local, ...) {
  n <- NROW(data)

  if (length(mean) != n && length(mean) != 1) {
    stop("mean vector must be length n or length 1 (recycled)")
  }

  if (spcov_params[["de"]] == 0 && (missing(randcov_params) || is.null(randcov_params))) {
    base_val <- replicate(samples, rnorm(n, sd = sqrt(spcov_params[["ie"]])))
  } else {
    ## convert sp to data frame (point geometry)
    attr_sp <- attr(class(data), "package")
    if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
      stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
    }

    ## convert sf to data frame (point geometry) (1d objects obsolete)
    ### see if data has sf class
    if (inherits(data, "sf")) {
      data <- suppressWarnings(sf::st_centroid(data))
      data <- sf_to_df(data)
      ### name xcoord ".xcoord" to be used later
      xcoord <- ".xcoord"
      ### name ycoord ".ycoord" to be used later
      ycoord <- ".ycoord"
    }

    # non standard evaluation for the x and y coordinates -- as.character()
    # right at capture normalizes both quoted ("x") and unquoted (x)
    # column-name references into a plain string
    xcoord <- as.character(substitute(xcoord))
    # replace null if necessary -- this missing() check runs before ycoord's
    # own capture below, so it still reflects the original argument correctly
    if (missing(ycoord)) {
      # 1-D data: fabricate a constant y-coordinate so the same 2-D distance
      # machinery below can be reused without a separate 1-D code path
      ycoord <- ".ycoord"
      data[[ycoord]] <- 0
    }
    ycoord <- as.character(substitute(ycoord))


    # storing x and y coordinate values
    xcoord_val <- data[[xcoord]]
    ycoord_val <- data[[ycoord]]

    # provide warning for this
    data$...response... <- seq(1, n)
    data$...xcoord... <- xcoord_val
    data$...ycoord... <- ycoord_val
    if ("extra" %in% names(spcov_params)) {
      spcov_init <- spcov_initial(
        spcov_type = class(spcov_params),
        de = spcov_params[["de"]],
        ie = spcov_params[["ie"]],
        range = spcov_params[["range"]],
        extra = spcov_params[["extra"]],
        rotate = spcov_params[["rotate"]],
        scale = spcov_params[["scale"]],
        known = "given"
      )
    } else {
      spcov_init <- spcov_initial(
        spcov_type = class(spcov_params),
        de = spcov_params[["de"]],
        ie = spcov_params[["ie"]],
        range = spcov_params[["range"]],
        rotate = spcov_params[["rotate"]],
        scale = spcov_params[["scale"]],
        known = "given"
      )
    }

    if (missing(randcov_params)) {
      randcov_params <- NULL
    } else {
      randcov_init <- randcov_initial(randcov_params, known = "given")
    }
    if (missing(partition_factor)) {
      partition_factor <- NULL
    }

    if (missing(local)) local <- NULL
    local_list <- get_local_list_simulation(local, n, data)

    if (local_list$approximation == "vecchia") {
      # vecchia: every location is simulated sequentially, conditional on
      # every earlier-simulated location so there is no base sample to
      # subset data down to at all (the object is built on all the data
      # and splm()'s own local = TRUE keeps that fit itself scalable; the
      # simulation is separately scalable via neighbor truncation)
      object <- splm(
        formula = ...response... ~ 1,
        data = data,
        spcov_initial = spcov_init,
        randcov_initial = randcov_init,
        partition_factor = partition_factor,
        xcoord = "...xcoord...",
        ycoord = "...ycoord...",
        local = TRUE
      )
      base_val <- get_sprnorm_vecchia(object, local_list, samples)
    } else {
      if (local_list$method_base != "all") {
        newdata <- lapply(local_list$index$new, function(x) data[x, , drop = FALSE])
        data <- data[local_list$index$base, , drop = FALSE]
        n <- NROW(data)
      }

      object <- splm(
        formula = ...response... ~ 1,
        data = data,
        spcov_initial = spcov_init,
        randcov_initial = randcov_init,
        partition_factor = partition_factor,
        xcoord = "...xcoord...",
        ycoord = "...ycoord...",
        local = TRUE
      )

      cov_lowchol_base <- t(chol(covmatrix(object)))
      base_val <- vapply(seq_len(samples), function(x) as.numeric(cov_lowchol_base %*% rnorm(n)), numeric(n))

      if (local_list$method_base != "all") {

        if (local_list$parallel) {
          cl <- parallel::makeCluster(local_list$ncores)
          new_val <- parLapply(cl, newdata, get_conditional_new_from_base, object, base_val, cov_lowchol_base, samples)
          cl <- parallel::stopCluster(cl)
        } else {
          new_val <- lapply(newdata, get_conditional_new_from_base, object, base_val, cov_lowchol_base, samples)
        }
        base_val <- rbind(base_val, do.call("rbind", new_val))
        index <- c(local_list$index$base, do.call("c", local_list$index$new))
        base_val <- base_val[order(index), , drop = FALSE]
      }
    }
  }

  base_val <- sweep(base_val, 1, mean, "+")


  if (samples == 1) {
    base_val <- as.vector(base_val)
  }
  base_val
}

#' @method sprnorm spherical
#' @export
sprnorm.spherical <- sprnorm.exponential

#' @method sprnorm gaussian
#' @export
sprnorm.gaussian <- sprnorm.exponential

#' @method sprnorm triangular
#' @export
sprnorm.triangular <- sprnorm.exponential

#' @method sprnorm circular
#' @export
sprnorm.circular <- sprnorm.exponential

#' @method sprnorm cubic
#' @export
sprnorm.cubic <- sprnorm.exponential

#' @method sprnorm pentaspherical
#' @export
sprnorm.pentaspherical <- sprnorm.exponential

#' @method sprnorm cosine
#' @export
sprnorm.cosine <- sprnorm.exponential

#' @method sprnorm wave
#' @export
sprnorm.wave <- sprnorm.exponential

#' @method sprnorm jbessel
#' @export
sprnorm.jbessel <- sprnorm.exponential

#' @method sprnorm gravity
#' @export
sprnorm.gravity <- sprnorm.exponential

#' @method sprnorm rquad
#' @export
sprnorm.rquad <- sprnorm.exponential

#' @method sprnorm magnetic
#' @export
sprnorm.magnetic <- sprnorm.exponential

#' @method sprnorm matern
#' @export
sprnorm.matern <- sprnorm.exponential

#' @method sprnorm cauchy
#' @export
sprnorm.cauchy <- sprnorm.exponential

#' @method sprnorm pexponential
#' @export
sprnorm.pexponential <- sprnorm.exponential

#' @rdname sprnorm
#' @method sprnorm none
#' @export
sprnorm.none <- function(spcov_params, mean = 0, samples = 1, data, randcov_params, partition_factor, ...) {
  n <- NROW(data)

  if (length(mean) != n && length(mean) != 1) {
    stop("mean vector must be length n or length 1 (recycled)")
  }

  dist_matrix <- diag(n)

  # compute the random effects covariance matrix
  if (missing(randcov_params)) {
    randcov_params <- NULL
    randcov_Zs <- NULL
  } else {
    names(randcov_params) <- get_randcov_names(reformulate(paste("(", names(randcov_params), ")", sep = "")))
    randcov_Zs <- get_randcov_Zs(data = data, names(randcov_params))
  }

  # partition matrix
  if (missing(partition_factor)) partition_factor <- NULL
  partition_matrix_val <- partition_matrix(partition_factor, data)

  # compute the covariance matrix
  cov_matrix_val <- cov_matrix(
    spcov_params, dist_matrix,
    randcov_params, randcov_Zs, partition_matrix_val
  )

  if (is.null(randcov_params)) {
    # with no spatial dependence and no random effects, the covariance matrix
    # is just ie * I, so drawing directly from rnorm() with sd = sqrt(ie) is
    # equivalent to (and much cheaper than) the general Cholesky route below
    sprnorm_val <- vapply(seq_len(samples), function(x) mean + rnorm(n, sd = sqrt(spcov_params[["ie"]])), numeric(n))
  } else {
    # transpose is lower triangular, needed for normal sim
    cov_matrix_lowchol <- t(chol(cov_matrix_val))
    # record sample sizes

    # simulate n random normal vectors
    sprnorm_val <- vapply(seq_len(samples), function(x) mean + as.numeric(cov_matrix_lowchol %*% rnorm(n)), numeric(n))
  }

  if (samples == 1) {
    sprnorm_val <- as.vector(sprnorm_val)
  }

  sprnorm_val
}


#' @rdname sprnorm
#' @method sprnorm ie
#' @export
sprnorm.ie <- sprnorm.none

#' @rdname sprnorm
#' @method sprnorm car
#' @export
sprnorm.car <- function(spcov_params, mean = 0, samples = 1, data, randcov_params, partition_factor, W, row_st = TRUE, M, ...) {
  n <- NROW(data)

  if (length(mean) != n && length(mean) != 1) {
    stop("mean vector must be length n or length 1 (recycled)")
  }


  # create distance matrix (if not provided) -- sf::st_intersects() assumes
  # units are nieghbors with themselves, so we need to set the diagonal of the
  # matrix equal to zero
  if (missing(W)) {
    ## convert sp to sf object
    attr_sp <- attr(class(data), "package")
    if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
      stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
    }
    W <- sf::st_intersects(data, sparse = FALSE)
    diag(W) <- 0
  }

  W <- 1 * Matrix::Matrix(W, sparse = TRUE)
  W_rowsums <- Matrix::rowSums(W)

  # make M if necessary
  if (row_st) {
    # under row standardization, M = diag(1 / rowSums(W)) is the matrix that
    # makes the CAR symmetry condition (I - range * W)^{-1} M symmetric hold
    if (!missing(M)) {
      warning("Overriding M when row_st = TRUE", call. = FALSE)
    }
    M <- 1 / W_rowsums # this has not been standardized
  } else {
    if (missing(M)) M <- rep(1, nrow(W)) # assume identity
  }

  if (row_st) {
    W_rowsums_val <- W_rowsums # make copy so rowsums are saved later
    # units with zero neighbors would otherwise divide by zero here; since
    # their entire row of W is already zero this substitution is a no-op
    W_rowsums_val[W_rowsums_val == 0] <- 1 # not a Matrix object so this subsetting is okay
    W <- W / W_rowsums_val
  }

  # verify the CAR symmetry condition holds so the resulting covariance
  # matrix (derived from (I - range * W)^{-1} M) is a valid, symmetric one
  if (inherits(spcov_params, "car") && !isSymmetric(as.matrix((Matrix(diag(nrow(W)), sparse = TRUE) - W) * 1 / M))) {
    stop("W and M must satisfy the CAR symmetry condition", call. = FALSE)
  }

  dist_matrix <- W

  # compute the random effects covariance matrix
  if (missing(randcov_params)) {
    randcov_params <- NULL
    randcov_Zs <- NULL
  } else {
    names(randcov_params) <- get_randcov_names(reformulate(paste("(", names(randcov_params), ")", sep = "")))
    randcov_Zs <- get_randcov_Zs(data = data, names(randcov_params))
  }

  # partition matrix
  if (missing(partition_factor)) partition_factor <- NULL
  partition_matrix_val <- partition_matrix(partition_factor, data)

  # compute the covariance matrix
  cov_matrix_val <- cov_matrix(
    spcov_params, dist_matrix,
    randcov_params, randcov_Zs, partition_matrix_val, M
  )

  # transpose is lower triangular, needed for normal sim
  cov_matrix_lowchol <- t(chol(cov_matrix_val))
  # record sample sizes

  # simulate n random normal vectors
  sprnorm_val <- vapply(seq_len(samples), function(x) mean + as.numeric(cov_matrix_lowchol %*% rnorm(n)), numeric(n))

  if (samples == 1) {
    sprnorm_val <- as.vector(sprnorm_val)
  }

  sprnorm_val
}

#' @method sprnorm sar
#' @export
sprnorm.sar <- sprnorm.car
