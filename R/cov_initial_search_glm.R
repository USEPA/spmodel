#' Search for initial covariance parameters
#'
#' @param spcov_initial_NA A spatial initial NA object
#' @param dispersion_initial_NA A dispersion initial NA object
#' @param ... Additional arguments passed to other methods
#' @noRd
cov_initial_search_glm <- function(spcov_initial_NA, dispersion_initial_NA, ...) {
  UseMethod("cov_initial_search_glm", spcov_initial_NA)
}
#' @export
cov_initial_search_glm.exponential <- function(spcov_initial_NA, dispersion_initial_NA, estmethod, data_object,
                                               dist_matrix_list, weights,
                                               randcov_initial_NA = NULL, esv_dotlist, ...) {
  # GLM counterpart to cov_initial_search.exponential(): builds the same de/ie/range
  # (and, with random effects, spatially-dominant/evenly-split/random-dominant)
  # starting-value grid, but also grids over a dispersion parameter and evaluates
  # the (Laplace-approximate) GLM objective via eval_grid_glm() instead of the
  # Gaussian objective, since sv-wls/sv-cl aren't applicable to non-Gaussian responses.
  # find ols sample variance
  s2 <- data_object$s2
  ns2 <- 1.2 * s2

  # find sets of starting values
  de <- c(0.1, 0.5, 0.9)
  # ie
  ie <- c(0.1, 0.5, 0.9)
  ## range
  range <- get_initial_range(class(spcov_initial_NA), data_object$max_halfdist) * c(0.5, 1.5)
  if (data_object$anisotropy) {
    ## rotate
    rotate <- c(0, 30 * pi / 180, 60 * pi / 180)
    ## scale
    scale <- c(0.25, 0.75, 1)
  } else {
    rotate <- 0
    scale <- 1
  }


  # find starting spatial grid
  spcov_grid <- build_de_ie_grid(de, ie, ns2, range = range, rotate = rotate, scale = scale)

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # replace with initial values
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }


  # take unique rows
  spcov_grid <- unique(spcov_grid)

  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # add dispersion: default candidate is 1 (i.e. no over/under-dispersion beyond
    # what the GLM family implies), overridden if the user fixed a value
    spcov_grid <- add_dispersion_column(spcov_grid, dispersion_initial_NA)


    # split
    min_params <- select_best_grid_point(spcov_grid, eval_grid_glm,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params

    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = NULL
    )
  } else {
    # randcov names
    randcov_names <- data_object$randcov_names
    cov_grid <- add_randcov_grids(spcov_grid, spcov_grid_init, ns2, spcov_initial_NA, randcov_initial_NA,
      nvar_spcov = 2
    )

    # add dispersion
    cov_grid <- add_dispersion_column(cov_grid, dispersion_initial_NA)

    min_params <- select_best_grid_point(cov_grid, eval_grid_glm,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the dispersion parameter
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = randcov_initial_NA
    )
  }

  # return the best parameters
  best_params
}

#' @export
cov_initial_search_glm.spherical <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.gaussian <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.triangular <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.circular <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.cubic <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.pentaspherical <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.cosine <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.wave <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.jbessel <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.gravity <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.rquad <- cov_initial_search_glm.exponential
#' @export
cov_initial_search_glm.magnetic <- cov_initial_search_glm.exponential

#' @export
cov_initial_search_glm.none <- function(spcov_initial_NA, dispersion_initial_NA, estmethod, data_object,
                                        dist_matrix_list, weights,
                                        randcov_initial_NA = NULL, esv_dotlist, ...) {
  # Unlike cov_initial_search.none() (the Gaussian version), this can't shortcut
  # straight to a closed-form ie estimate even with no random effects, because the
  # dispersion parameter still needs to be grid-searched via eval_grid_glm()
  # find ols sample variance
  s2 <- data_object$s2
  ns2 <- 1.2 * s2

  # find sets of starting values
  ## de
  de <- 0
  ## ie
  ie <- 1
  ## range
  range <- get_initial_range(class(spcov_initial_NA), NULL)
  ## rotate
  rotate <- 0
  ## scale
  scale <- 1


  # find starting spatial grid
  spcov_grid <- build_de_ie_grid(de, ie, ns2, range = range, rotate = rotate, scale = scale)

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # replace with initial values
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }

  # take unique rows
  spcov_grid <- unique(spcov_grid)

  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # add dispersion
    spcov_grid <- add_dispersion_column(spcov_grid, dispersion_initial_NA)

    # split
    min_params <- select_best_grid_point(spcov_grid, eval_grid_glm,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params

    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = NULL
    )
  } else {
    # randcov vars names
    randcov_names <- data_object$randcov_names
    # none/ie has no meaningful de/ie spread to filter the evenly-dominated or
    # random-dominated grids on (de is always 0), so every row is kept
    all_rows <- rep(TRUE, NROW(spcov_grid_init))
    cov_grid <- add_randcov_grids(spcov_grid, spcov_grid_init, ns2, spcov_initial_NA, randcov_initial_NA,
      nvar_spcov = 1, evencov_filter = all_rows, randcov_filter = all_rows
    )

    # add dispersion
    cov_grid <- add_dispersion_column(cov_grid, dispersion_initial_NA)

    min_params <- select_best_grid_point(cov_grid, eval_grid_glm,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the dispersion parameter
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = randcov_initial_NA
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search_glm.ie <- cov_initial_search_glm.none

#' @export
cov_initial_search_glm.matern <- function(spcov_initial_NA, dispersion_initial_NA, estmethod, data_object,
                                          dist_matrix_list, weights,
                                          randcov_initial_NA = NULL, esv_dotlist, ...) {
  # Same grid-search strategy as cov_initial_search_glm.exponential() (see comments
  # there), extended with a starting grid for the extra shape/smoothness parameter
  # shared by matern-family correlation functions (matern, cauchy, pexponential)
  # find ols sample variance
  s2 <- data_object$s2
  ns2 <- 1.2 * s2


  # find sets of starting values
  de <- c(0.1, 0.5, 0.9)
  # ie
  ie <- c(0.1, 0.5, 0.9)
  ## range
  range <- get_initial_range(class(spcov_initial_NA), data_object$max_halfdist) * c(0.5, 1.5)
  ## extra
  extra <- get_initial_extra(class(spcov_initial_NA)) * c(0.5, 2)
  if (data_object$anisotropy) {
    ## rotate
    rotate <- c(0, 30 * pi / 180, 60 * pi / 180)
    ## scale
    scale <- c(0.25, 0.75, 1)
  } else {
    rotate <- 0
    scale <- 1
  }


  # find starting spatial grid
  spcov_grid <- build_de_ie_grid(de, ie, ns2, range = range, extra = extra, rotate = rotate, scale = scale)

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # replace with initial values
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }

  # take unique rows
  spcov_grid <- unique(spcov_grid)

  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # add dispersion
    spcov_grid <- add_dispersion_column(spcov_grid, dispersion_initial_NA)

    # split
    min_params <- select_best_grid_point(spcov_grid, eval_grid_glm,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    spcov_params <- min_params[c("de", "ie", "range", "extra", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params

    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = NULL
    )
  } else {
    # randcov names
    randcov_names <- data_object$randcov_names
    cov_grid <- add_randcov_grids(spcov_grid, spcov_grid_init, ns2, spcov_initial_NA, randcov_initial_NA,
      nvar_spcov = 2,
      randcov_filter = spcov_grid_init$de == spcov_grid_init$ie &
        spcov_grid_init$range == min(spcov_grid_init$range) &
        spcov_grid_init$extra == min(spcov_grid_init$extra)
    )
    # add dispersion
    cov_grid <- add_dispersion_column(cov_grid, dispersion_initial_NA)

    min_params <- select_best_grid_point(cov_grid, eval_grid_glm,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = dist_matrix_list
    )
    spcov_params <- min_params[c("de", "ie", "range", "extra", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the dispersion parameter
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = randcov_initial_NA
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search_glm.cauchy <- cov_initial_search_glm.matern
#' @export
cov_initial_search_glm.pexponential <- cov_initial_search_glm.matern

#' @export
cov_initial_search_glm.car <- function(spcov_initial_NA, dispersion_initial_NA, estmethod, data_object,
                                       dist_matrix_list, randcov_initial_NA = NULL, ...) {
  # GLM counterpart to cov_initial_search.car(): "range" candidates are spread
  # across the valid autocorrelation (rho) range rather than spatial distances, and
  # a dispersion grid is added (see cov_initial_search_glm.exponential() comments)
  # find ols sample variance
  s2 <- data_object$s2
  ns2 <- 1.2 * s2

  # store W as dist_matrix
  # MAKE IT CLEAR DIST_MATRIX_LIST IS NOT A LIST
  W <- dist_matrix_list

  # find sets of starting values
  ## de
  de <- c(0.1, 0.5, 0.9)
  ## ie
  ie <- c(0.1, 0.5, 0.9)
  ## range
  rho_length <- data_object$rho_ub - data_object$rho_lb
  range <- c(
    data_object$rho_lb + 0.01 * rho_length,
    mean(c(data_object$rho_lb, data_object$rho_ub)),
    data_object$rho_ub - 0.01 * rho_length
  )
  # find starting spatial grid
  spcov_grid <- build_de_ie_grid(de, ie, ns2, range = range)
  spcov_grid$extra <- spcov_grid$de

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # replace with initial values
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }

  # take unique rows
  spcov_grid <- unique(spcov_grid)

  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # add dispersion
    spcov_grid <- add_dispersion_column(spcov_grid, dispersion_initial_NA)

    min_params <- select_best_grid_point(spcov_grid, eval_grid_glm,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = W
    )
    spcov_params <- min_params[c("de", "ie", "range", "extra")]
    spcov_initial_NA$initial <- spcov_params
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = NULL
    )
  } else {
    # randcov vars names
    randcov_names <- data_object$randcov_names
    cov_grid <- add_randcov_grids(spcov_grid, spcov_grid_init, ns2, spcov_initial_NA, randcov_initial_NA,
      nvar_spcov = 2, scale_cols = c("de", "ie", "extra")
    )

    # add dispersion
    cov_grid <- add_dispersion_column(cov_grid, dispersion_initial_NA)

    min_params <- select_best_grid_point(cov_grid, eval_grid_glm,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      family = class(dispersion_initial_NA), estmethod = estmethod, dist_matrix_list = W
    )
    spcov_params <- min_params[c("de", "ie", "range", "extra")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the dispersion parameter
    dispersion_params <- min_params[c("dispersion")]
    dispersion_initial_NA$initial <- dispersion_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA,
      dispersion_initial_val = dispersion_initial_NA, randcov_initial_val = randcov_initial_NA
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search_glm.sar <- cov_initial_search_glm.car

#' Evaluate the (non-spatial-random-effect) grid-search objective for one grid point
#'
#' @param cov_grid_split A one-row split of the covariance parameter grid
#' @param data_object The data object
#' @param spcov_type The spatial covariance type
#' @param family The GLM family
#' @param estmethod The estimation method
#' @param dist_matrix_list A list of distance matrices
#'
#' @return The (REML/ML/sv-wls/sv-cl, as determined by \code{estmethod}) objective
#'   function value at this grid point, used by \code{cov_initial_search_glm()}
#'   methods to choose starting values
#'
#' @noRd
eval_grid_glm <- function(cov_grid_split, data_object, spcov_type,
                          family, estmethod, dist_matrix_list) {
  # convert list structure to a vector
  cov_grid <- unlist(cov_grid_split)

  # find spatial covariance parameter vector
  spcov_grid <- cov_grid[c("de", "ie", "range", "extra", "rotate", "scale")]
  spcov_grid <- spcov_grid[!is.na(spcov_grid)]
  spcov_params_val <- do.call("spcov_params", c(list(spcov_type = spcov_type), as.list(spcov_grid)))

  # dispersion parameter vector
  dispersion_params_val <- dispersion_params(family = data_object$family, dispersion = cov_grid[["dispersion"]])

  # find REML or ML objective function value (via the Laplace approximation, since
  # GLM responses generally have no closed-form marginal likelihood)

  # incorporate random effects if necessary
  if (is.null(data_object$randcov_initial)) {
    randcov_params_val <- NULL
  } else {
    randcov_names <- data_object$randcov_names
    randcov_params_val <- randcov_params(cov_grid[randcov_names], nm = randcov_names)
  }

  # incorporate anisotropy if necessary: as in eval_grid(), rotate is only
  # identifiable modulo pi, so both candidate angles are evaluated and the smaller
  # (better) objective value is kept
  if (data_object$anisotropy) {
    objval <- resolve_rotation_ambiguity(spcov_params_val, randcov_params_val, data_object, estmethod,
      dispersion_params_val = dispersion_params_val
    )
  } else {
    # compute relevant products
    lapll_prods <- laploglik_products(
      spcov_params_val, dispersion_params_val, data_object, estmethod,
      dist_matrix_list, randcov_params_val
    )

    # find -2loglik
    objval <- get_minustwolaploglik(lapll_prods, estmethod, data_object$n,
      data_object$p,
      spcov_profiled = FALSE, randcov_profiled = FALSE
    )
  }
  objval
}
