#' Create a Spatial Decorrelation Transformation Grid
#'
#' @description
#'  Create a spatial decorrelation transformation grid of initial parameters to be
#'   evaluated via a grid search.
#'
#' @inheritParams decorrelate
#'
#' @param dense_grid A logical
#'   which controls the density of the constructed grid to be evaluated. If
#'   \code{dense_grid} is \code{TRUE}, a denser grid is used. If \code{dense_grid}
#'   is \code{FALSE}, a sparser grid is used. By default, \code{dense_grid}
#'   is \code{FALSE} when the sample size is greater than 5,000 and \code{TRUE}
#'   otherwise.
#'
#' @return A grid of spatial decorrelation parameters stored as a \code{data.frame}.
#'
#' @export
#'
#' @examples
#' decorrelate_grid(log_cond ~ temp, data = lake, spcov_type = "exponential")
decorrelate_grid <- function(formula, data, spcov_type, spcov_params, xcoord, ycoord, anisotropy = FALSE, random, randcov_params, dense_grid) {

  # non standard evaluation for x and y coordinates (only meaningful at this,
  # the direct calling frame -- decorrelate_grid_internal() receives the
  # already-substituted value and must not re-substitute)
  xcoord <- substitute(xcoord)
  ycoord <- substitute(ycoord)

  if (missing(spcov_params)) spcov_params <- NULL
  if (!is.null(spcov_params)) spcov_type <- class(spcov_params)
  if (missing(random)) random <- NULL
  if (missing(randcov_params)) randcov_params <- NULL

  if (missing(dense_grid)) {
    if (NROW(data) <= 5000) {
      dense_grid <- TRUE
    } else {
      dense_grid <- FALSE
    }
  }

  decorrelate_grid_internal(
    formula = formula,
    data = data,
    spcov_type = spcov_type,
    spcov_params = spcov_params,
    xcoord = xcoord,
    ycoord = ycoord,
    anisotropy = anisotropy,
    random = random,
    randcov_params = randcov_params,
    dense_grid = dense_grid,
    add_iid = TRUE,
    warn = TRUE
  )
}

# Shared worker behind decorrelate_grid() (called with add_iid = TRUE, warn = TRUE)
# and decorrelate_initial_search() (called once per grid-search training split,
# typically with warn = FALSE to avoid repeating the same geometry-coercion
# warning across replications/folds, and add_iid depending on whether an "iid"
# baseline comparison row is already implied by the caller's inputs).
#
# xcoord/ycoord must already be resolved (not NSE symbols to substitute) and
# dense_grid must already be resolved -- both are the caller's responsibility,
# since substitute() only works meaningfully in decorrelate_grid()'s own frame.
decorrelate_grid_internal <- function(formula, data, spcov_type, spcov_params, xcoord, ycoord, anisotropy = FALSE, random, randcov_params, dense_grid, add_iid, warn) {



  if (missing(spcov_params)) spcov_params <- NULL
  if (!is.null(spcov_params)) spcov_type <- class(spcov_params)
  if (missing(random)) random <- NULL
  if (missing(randcov_params)) randcov_params <- NULL
  if (!is.null(randcov_params)) {
    # overwrite random if randcov_params provided
    random <- reformulate(names(randcov_params))
  }

  # find ols sample variance
  lmod <- lm(formula, data)
  s2 <- summary(lmod)$sigma^2
  ns2 <- 1.2 * s2

  # find sets of starting values
  ## de
  # de <- c(0.1, 0.5, 0.9)
  if (dense_grid) {
    de <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  } else {
    de <- c(0.5, 0.95)
  }

  ## ie
  # ie <- c(0.1, 0.5, 0.9)
  if (dense_grid) {
    ie <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  } else {
    ie <- c(0.05, 0.5)
  }
  ## range
  if (inherits(data, "sf")) {
    if (!spcov_type %in% c("none", "ie") && any(sf::st_geometry_type(data) != "POINT")) {
      if (warn) {
        warning("At least one geometry type in data is not equal to \"POINT\". Attempting to coerce all non-\"POINT\" geometries to \"POINT\" geometries via their centroids using sf::st_centroid().", call. = FALSE)
      }
    }
    data_sf <- suppressWarnings(sf::st_centroid(data))
    # store as data frame
    data <- sf_to_df(data_sf)
    ## name xcoord ".xcoord" to be used later
    xcoord <- ".xcoord"
    ## name ycoord ".ycoord" to be used later
    ycoord <- ".ycoord"
  }
  # storing max halfdist
  x_range <- range(data[[xcoord]])
  if (spcov_type %in% c("triangular", "cosine")) {
    data[[ycoord]] <- 0
  }
  y_range <- range(data[[ycoord]])
  max_halfdist <- sqrt((max(x_range) - min(x_range))^2 + (max(y_range) - min(y_range))^2) / 2
  range <- get_initial_range(spcov_type, max_halfdist) * c(0.5, 1.5)
  ## anisotropy
  if (anisotropy) {
    ## rotate
    # rotate <- c(0, 30 * pi / 180, 60 * pi / 180)
    if (dense_grid) {
      rotate <- c(0, 45, 90, 135) * pi / 180
    } else {
      rotate <- c(0, 90) * pi / 180
    }

    ## scale
    # scale <- c(0.25, 0.75, 1)
    if (dense_grid) {
      scale <- c(0.25, 0.5, 0.75, 1)
    } else {
      scale <- c(0.5, 1)
    }
  } else {
    ## rotate
    rotate <- 0
    ## scale
    scale <- 1
  }

  # find starting spatial grid
  spcov_grid <- expand.grid(de = de, ie = ie, range = range, rotate = rotate, scale = scale)
  if (spcov_type %in% c("matern", "cauchy", "pexponential")) {
    ## extra (if applicable)
    extra <- get_initial_extra(spcov_type) * c(0.5, 2)
    # recalculating less efficient but more readable
    spcov_grid <- expand.grid(de = de, ie = ie, range = range, extra = extra, rotate = rotate, scale = scale)
  }
  spcov_grid <- spcov_grid[spcov_grid$de + spcov_grid$ie == 1, , drop = FALSE]
  spcov_grid[, c("de", "ie")] <- ns2 * spcov_grid[, c("de", "ie")]

  # anisotropy correction
  spcov_grid$rotate[spcov_grid$scale == 1] <- 0

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid
  # take unique rows
  spcov_grid <- unique(spcov_grid)

  if (!missing(random) && !is.null(random)) {
    randcov_names <- get_randcov_names(random)
    # find number of random effects
    nvar_randcov <- length(randcov_names)

    # spatially dominant grid
    ## spatial components 90% of variance
    spcov_grid[, c("de", "ie")] <- 0.9 * spcov_grid[, c("de", "ie")]
    for (x in randcov_names) {
      spcov_grid[, x] <- 0.1 * ns2 / nvar_randcov
    }

    # Evenly dominated grid
    ## find number of spatial variance parameters
    nvar_spcov <- 2
    ## find all variance parameters
    nvar_cov <- nvar_spcov + nvar_randcov
    ## keep only spatial cases with variance spread out
    evencov_grid <- spcov_grid_init[spcov_grid_init$de == spcov_grid_init$ie, , drop = FALSE]
    ## scale to incorporate overall variance
    evencov_grid[, c("de", "ie")] <- nvar_spcov / nvar_cov * evencov_grid[, c("de", "ie")]
    # find unique combinations
    evencov_grid <- unique(evencov_grid)
    # add random effects
    for (x in randcov_names) {
      evencov_grid[, x] <- 1 / nvar_cov * ns2
    }

    # random dominated grid
    ## spatial component 10%
    ## include all range parameter values for decorrelate grid search
    randcov_grid_spcov <- spcov_grid_init[spcov_grid_init$de == spcov_grid_init$ie & spcov_grid_init$range == min(spcov_grid_init$range), , drop = FALSE]
    randcov_grid_spcov[, c("de", "ie")] <- 0.1 * randcov_grid_spcov[, c("de", "ie")]

    # random dominant grid
    randcov_grid_randcov <- as.data.frame(as.list(rep(1 / nvar_randcov, nvar_randcov)))
    names(randcov_grid_randcov) <- randcov_names
    ## if there is more than one random effect, split it up into relevant scenarios
    if (nvar_randcov > 1) {
      ## set 0.1 for all proportions
      extra_grid_randcov <- lapply(seq_len(nvar_randcov), function(x) rep(0.1 * 1 / (nvar_randcov - 1), nvar_randcov))
      extra_grid_randcov <- do.call("rbind", extra_grid_randcov)
      ## fill in 0.9 for one proportion in each row
      diag(extra_grid_randcov) <- 0.9
      extra_grid_randcov <- as.data.frame(extra_grid_randcov)
      names(extra_grid_randcov) <- randcov_names
    } else {
      extra_grid_randcov <- NULL
    }
    ## give the random effects 90% of the variance
    randcov_grid_randcov <- 0.9 * ns2 * rbind(randcov_grid_randcov, extra_grid_randcov)

    ## bind together and replicate
    randcov_grid_spcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_randcov), randcov_grid_spcov, simplify = FALSE))
    randcov_grid_randcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_spcov), randcov_grid_randcov, simplify = FALSE))
    randcov_grid <- cbind(randcov_grid_spcov_rep, randcov_grid_randcov_rep)

    # bind together all grids
    cov_grid <- rbind(spcov_grid, evencov_grid, randcov_grid)
  } else {
    cov_grid <- spcov_grid
  }
  cov_grid$spcov_type <- spcov_type
  ncols <- NCOL(cov_grid)
  cov_grid <- cov_grid[, c(ncols, seq(1, ncols - 1))]
  if (spcov_type %in% c("none", "ie")) {
    cov_grid$de <- 0
    if (spcov_type == "none") {
      cov_grid$ie <- ns2
      # cov_grid$ie <- 1
    } else if (spcov_type == "ie") {
      cov_grid$ie <- ns2
    }
    cov_grid$range <- Inf
    cov_grid$rotate <- 0
    cov_grid$scale <- 1
    anisotropy <- FALSE
  }
  cov_grid$spcov_type[cov_grid$de == 0] <- "none"
  # if (!anisotropy) {
  #   remove_cols <- which(names(cov_grid) %in% c("rotate", "scale"))
  #   cov_grid <- cov_grid[, -remove_cols, drop = FALSE]
  # }


  if (!is.null(spcov_params)) {
    for (x in names(spcov_params)) {
      cov_grid[, x] <- spcov_params[[x]]
    }
  }
  if (!is.null(randcov_params)) {
    randcov_names <- get_randcov_names(random)
    names(randcov_params) <- randcov_names
    for (x in randcov_names) {
      cov_grid[, x] <- randcov_params[[x]]
    }
  }

  # create iid grid
  if (add_iid) {
    iid_grid <- as.data.frame(matrix(0, nrow = 1, ncol = NCOL(cov_grid)))
    names(iid_grid) <- names(cov_grid)
    iid_grid$spcov_type <- "none"
    iid_grid$ie <- 1
    iid_grid$range <- Inf
    iid_grid$scale <- 1
    cov_grid <- rbind(cov_grid, iid_grid)
  }


  # return cov grid
  cov_grid <- unique(cov_grid)
  row.names(cov_grid) <- as.character(seq(1, NROW(cov_grid)))
  cov_grid
}
