#' Search for initial covariance parameters
#'
#' @param spcov_initial_NA A spatial initial NA object
#' @param ... Additional arguments passed to other methods
#' @noRd
cov_initial_search <- function(spcov_initial_NA, ...) {
  UseMethod("cov_initial_search", spcov_initial_NA)
}
#' @export
cov_initial_search.exponential <- function(spcov_initial_NA, estmethod, data_object,
                                           dist_matrix_list, weights,
                                           randcov_initial_NA = NULL, esv_dotlist, ...) {
  # Rather than starting REML/ML/sv-wls/sv-cl optimization from one arbitrary
  # starting point (which risks a poor local optimum), this builds a small grid of
  # plausible starting values -- combinations of "de" (dependent/partial-sill
  # variance, the spatially structured variance) and "ie" (independent-error/nugget
  # variance) that sum to a fixed variance budget, crossed with candidate ranges
  # (and rotate/scale under anisotropy) -- evaluates the objective at each grid
  # point (see eval_grid()), and returns the best one as the optimizer's start.
  # find ols sample variance
  s2 <- data_object$s2
  # inflate the OLS variance slightly (by 20%) as a rough total-variance budget to
  # split across the de/ie (and, further below, random effect) components, since
  # OLS residual variance tends to understate the true total variance once spatial
  # dependence is accounted for
  ns2 <- 1.2 * s2

  # find sets of starting values
  # ## de
  # de <- ns2 * c(0.1, 0.5, 0.9)
  # ## ie
  # ie <- ns2 * c(0.1, 0.5, 0.9)
  # de: proportions of ns2 attributed to spatially structured variance
  de <- c(0.1, 0.5, 0.9)
  # ie: proportions of ns2 attributed to the nugget/independent error
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


  # find starting spatial grid (keeping only combinations where de + ie
  # proportions sum to 1, i.e. de and ie together exhaust the full variance
  # budget with no double counting)
  spcov_grid <- build_de_ie_grid(de, ie, ns2, range = range, rotate = rotate, scale = scale)

  # save initial state (used with random effects)
  spcov_grid_init <- spcov_grid

  # any parameter the user fixed (not NA in spcov_initial_NA) overrides the grid
  # value in every row, since there's no point searching over a value that's fixed
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }

  # take unique rows (fixing parameters above can collapse several grid rows to duplicates)
  spcov_grid <- unique(spcov_grid)
  # for sv-wls, the objective (evaluated per grid point in eval_grid()) is a
  # weighted least-squares fit to this empirical semivariogram, so it only needs to
  # be computed once, up front, rather than inside the grid loop
  sv_wls_result <- precompute_sv_wls(estmethod, data_object, spcov_initial_NA, dist_matrix_list, esv_dotlist)
  esv_val <- sv_wls_result$esv_val
  dist_matrix_list <- sv_wls_result$dist_matrix_list

  # for sv-cl (composite likelihood via pairwise squared differences, Curriero &
  # Lele 1999), precompute the vector of pairwise distances and squared OLS-residual
  # differences once, up front, since the composite-likelihood objective in
  # eval_grid() is a function of these fixed vectors, not of the grid point itself
  sv_cl_result <- precompute_sv_cl(estmethod, data_object, spcov_initial_NA, dist_matrix_list)
  dist_vector <- sv_cl_result$dist_vector
  residual_vector2 <- sv_cl_result$residual_vector2
  dist_matrix_list <- sv_cl_result$dist_matrix_list


  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # evaluate the objective function at every grid point, keeping the
    # combination that minimizes it as the optimizer's starting point
    min_params <- select_best_grid_point(spcov_grid, eval_grid,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      estmethod = estmethod, dist_matrix_list = dist_matrix_list,
      weights = weights, esv = esv_val,
      dist_vector = dist_vector,
      residual_vector2 = residual_vector2
    )
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params

    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA, randcov_initial_val = NULL, esv = esv_val,
      dist_vector = dist_vector, residual_vector2 = residual_vector2
    )
  } else {
    # With random effects present, the total variance can plausibly be split many
    # ways between spatial dependence, nugget, and each random effect. Rather than
    # crossing every combination (a computational burden), three
    # representative grids are built below -- one where spatial variance dominates,
    # one where variance is spread evenly, and one where the random effect(s)
    # dominate -- and unioned together as the candidate starting points.
    randcov_names <- data_object$randcov_names
    cov_grid <- add_randcov_grids(spcov_grid, spcov_grid_init, ns2, spcov_initial_NA, randcov_initial_NA,
      nvar_spcov = 2
    )

    min_params <- select_best_grid_point(cov_grid, eval_grid,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      estmethod = estmethod, dist_matrix_list = dist_matrix_list,
      weights = weights, esv = esv_val,
      dist_vector = dist_vector,
      residual_vector2 = residual_vector2
    )
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA, randcov_initial_val = randcov_initial_NA, esv = esv_val,
      dist_vector = dist_vector, residual_vector2 = residual_vector2
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search.spherical <- cov_initial_search.exponential
#' @export
cov_initial_search.gaussian <- cov_initial_search.exponential
#' @export
cov_initial_search.triangular <- cov_initial_search.exponential
#' @export
cov_initial_search.circular <- cov_initial_search.exponential
#' @export
cov_initial_search.cubic <- cov_initial_search.exponential
#' @export
cov_initial_search.pentaspherical <- cov_initial_search.exponential
#' @export
cov_initial_search.cosine <- cov_initial_search.exponential
#' @export
cov_initial_search.wave <- cov_initial_search.exponential
#' @export
cov_initial_search.jbessel <- cov_initial_search.exponential
#' @export
cov_initial_search.gravity <- cov_initial_search.exponential
#' @export
cov_initial_search.rquad <- cov_initial_search.exponential
#' @export
cov_initial_search.magnetic <- cov_initial_search.exponential

#' @export
cov_initial_search.none <- function(spcov_initial_NA, estmethod, data_object,
                                    dist_matrix_list, weights,
                                    randcov_initial_NA = NULL, esv_dotlist, ...) {
  # "none" (and "ie", aliased below) means no spatial covariance structure at all --
  # the only variance component is the nugget/independent error (plus any random
  # effects). With no random effects there's nothing to grid-search over: the OLS
  # sample variance is directly the ie estimate, so this returns immediately.
  # find ols sample variance
  s2 <- data_object$s2

  # exit if no random effects
  if (is.null(randcov_initial_NA)) {
    spcov_initial_NA$initial["ie"] <- s2
    best_params <- list(
      spcov_initial_val = spcov_initial_NA, randcov_initial_val = NULL, esv = NULL,
      dist_vector = NULL, residual_vector2 = NULL
    )
    return(best_params)
  }

  # do other stuff
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

  # compute empirical semivariogram
  sv_wls_result <- precompute_sv_wls(estmethod, data_object, spcov_initial_NA, dist_matrix_list, esv_dotlist)
  esv_val <- sv_wls_result$esv_val
  dist_matrix_list <- sv_wls_result$dist_matrix_list

  # find relevant quantities for composite likelihood
  sv_cl_result <- precompute_sv_cl(estmethod, data_object, spcov_initial_NA, dist_matrix_list)
  dist_vector <- sv_cl_result$dist_vector
  residual_vector2 <- sv_cl_result$residual_vector2
  dist_matrix_list <- sv_cl_result$dist_matrix_list


  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # split
    min_params <- select_best_grid_point(spcov_grid, eval_grid,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      estmethod = estmethod, dist_matrix_list = dist_matrix_list,
      weights = weights, esv = esv_val,
      dist_vector = dist_vector,
      residual_vector2 = residual_vector2
    )
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA, randcov_initial_val = NULL, esv = esv_val,
      dist_vector = dist_vector, residual_vector2 = residual_vector2
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

    min_params <- select_best_grid_point(cov_grid, eval_grid,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      estmethod = estmethod, dist_matrix_list = dist_matrix_list,
      weights = weights, esv = esv_val,
      dist_vector = dist_vector,
      residual_vector2 = residual_vector2
    )
    spcov_params <- min_params[c("de", "ie", "range", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA, randcov_initial_val = randcov_initial_NA, esv = esv_val,
      dist_vector = dist_vector, residual_vector2 = residual_vector2
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search.ie <- cov_initial_search.none

#' @export
cov_initial_search.matern <- function(spcov_initial_NA, estmethod, data_object,
                                      dist_matrix_list, weights,
                                      randcov_initial_NA = NULL, esv_dotlist, ...) {
  # Same grid-search strategy as cov_initial_search.exponential() (see comments
  # there), extended with a starting grid for the extra shape/smoothness parameter
  # that matern-family correlation functions (matern, cauchy, pexponential) have in
  # addition to de/ie/range/rotate/scale
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

  # compute empirical semivariogram
  sv_wls_result <- precompute_sv_wls(estmethod, data_object, spcov_initial_NA, dist_matrix_list, esv_dotlist)
  esv_val <- sv_wls_result$esv_val
  dist_matrix_list <- sv_wls_result$dist_matrix_list

  # find relevant quantities for composite likelihood
  sv_cl_result <- precompute_sv_cl(estmethod, data_object, spcov_initial_NA, dist_matrix_list)
  dist_vector <- sv_cl_result$dist_vector
  residual_vector2 <- sv_cl_result$residual_vector2
  dist_matrix_list <- sv_cl_result$dist_matrix_list


  # perform search if no random effects
  if (is.null(randcov_initial_NA)) {
    # split
    min_params <- select_best_grid_point(spcov_grid, eval_grid,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      estmethod = estmethod, dist_matrix_list = dist_matrix_list,
      weights = weights, esv = esv_val,
      dist_vector = dist_vector,
      residual_vector2 = residual_vector2
    )
    spcov_params <- min_params[c("de", "ie", "range", "extra", "rotate", "scale")]
    spcov_initial_NA$initial <- spcov_params
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA, randcov_initial_val = NULL, esv = esv_val,
      dist_vector = dist_vector, residual_vector2 = residual_vector2
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

    min_params <- select_best_grid_point(cov_grid, eval_grid,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      estmethod = estmethod, dist_matrix_list = dist_matrix_list,
      weights = weights, esv = esv_val,
      dist_vector = dist_vector,
      residual_vector2 = residual_vector2
    )
    spcov_params <- min_params[c("de", "ie", "range", "extra", "rotate", "scale")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(
      spcov_initial_val = spcov_initial_NA, randcov_initial_val = randcov_initial_NA, esv = esv_val,
      dist_vector = dist_vector, residual_vector2 = residual_vector2
    )
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search.cauchy <- cov_initial_search.matern
#' @export
cov_initial_search.pexponential <- cov_initial_search.matern

#' @export
cov_initial_search.car <- function(spcov_initial_NA, estmethod, data_object,
                                   dist_matrix_list, randcov_initial_NA = NULL, ...) {
  # Areal (CAR, and SAR aliased below) autoregressive analogue of the geostatistical
  # grid search above: "range" here is really the autocorrelation parameter rho, so
  # its candidate values are spread across rho's valid range (rho_lb, rho_ub) --
  # the bounds within which the CAR/SAR precision matrix stays positive definite --
  # rather than across spatial distances, and there is no empirical-semivariogram
  # (sv-wls/sv-cl) branch since those aren't defined for areal data.
  # find ols sample variance
  # obdata <- data_object$data[data_object$observed_index, , drop = FALSE]
  # s2 <- summary(lm(data_object$formula, obdata))$sigma^2
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
  ## range: candidate autocorrelation (rho) values, kept strictly inside (rho_lb, rho_ub)
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
    # split
    min_params <- select_best_grid_point(spcov_grid, eval_grid,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      estmethod = estmethod, dist_matrix_list = W
    )
    spcov_params <- min_params[c("de", "ie", "range", "extra")]
    spcov_initial_NA$initial <- spcov_params
    # return the best parameters
    best_params <- list(spcov_initial_val = spcov_initial_NA, randcov_initial_val = NULL)
  } else {
    # randcov vars names
    randcov_names <- data_object$randcov_names
    cov_grid <- add_randcov_grids(spcov_grid, spcov_grid_init, ns2, spcov_initial_NA, randcov_initial_NA,
      nvar_spcov = 2, scale_cols = c("de", "ie", "extra")
    )

    min_params <- select_best_grid_point(cov_grid, eval_grid,
      data_object = data_object, spcov_type = class(spcov_initial_NA),
      estmethod = estmethod, dist_matrix_list = W
    )
    spcov_params <- min_params[c("de", "ie", "range", "extra")]
    # return the spatial parameters
    spcov_initial_NA$initial <- spcov_params
    # return the random effect parameters
    randcov_initial_NA$initial <- randcov_params(min_params[randcov_names])
    # return the best parameters
    best_params <- list(spcov_initial_val = spcov_initial_NA, randcov_initial_val = randcov_initial_NA)
  }
  # return the best parameters
  best_params
}

#' @export
cov_initial_search.sar <- cov_initial_search.car

#' Evaluate the grid-search objective for one grid point
#'
#' @param cov_grid_split A one-row split of the covariance parameter grid
#' @param data_object The data object
#' @param spcov_type The spatial covariance type
#' @param estmethod The estimation method
#' @param dist_matrix_list A list of distance matrices
#' @param weights Semivariogram weights (used when \code{estmethod} is \code{"sv-wls"})
#' @param esv An empirical semivariogram object (used when \code{estmethod} is \code{"sv-wls"})
#' @param dist_vector A distance vector (used when \code{estmethod} is \code{"sv-cl"})
#' @param residual_vector2 A vector of squared residuals (used when \code{estmethod} is \code{"sv-cl"})
#'
#' @return The (REML/ML/sv-wls/sv-cl, as determined by \code{estmethod}) objective
#'   function value at this grid point, used by \code{cov_initial_search()}
#'   methods to choose starting values
#'
#' @noRd
eval_grid <- function(cov_grid_split, data_object, spcov_type,
                      estmethod, dist_matrix_list,
                      weights, esv, dist_vector, residual_vector2) {
  # convert list structure to a vector
  cov_grid <- unlist(cov_grid_split)

  # find spatial covariance parameter vector
  spcov_grid <- cov_grid[c("de", "ie", "range", "extra", "rotate", "scale")]
  spcov_grid <- spcov_grid[!is.na(spcov_grid)]
  spcov_params_val <- do.call("spcov_params", c(list(spcov_type = spcov_type), as.list(spcov_grid)))

  # dispatch to the objective matching the requested estimation method: -2 times
  # the (restricted) Gaussian log-likelihood for reml/ml, or the sv-wls/sv-cl loss
  # otherwise -- lower is better in every case, so the caller just takes which.min()
  # find REML or ML objective function value
  if (estmethod %in% c("reml", "ml")) {
    # incorporate random effects if necessary
    if (is.null(data_object$randcov_initial)) {
      randcov_params_val <- NULL
    } else {
      randcov_names <- data_object$randcov_names
      randcov_params_val <- randcov_params(cov_grid[randcov_names], nm = randcov_names)
    }

    # incorporate anisotropy if necessary: the likelihood is evaluated at both
    # candidate angles (rotate and abs(pi - rotate))
    # and the better one is kept
    if (data_object$anisotropy) {
      objval <- resolve_rotation_ambiguity(spcov_params_val, randcov_params_val, data_object, estmethod)
    } else {
      # compute relevant products
      gll_prods <- gloglik_products(
        spcov_params_val, data_object, estmethod,
        dist_matrix_list, randcov_params_val
      )

      # find -2loglik
      objval <- get_minustwologlik(gll_prods, estmethod, data_object$n,
        data_object$p,
        spcov_profiled = FALSE, randcov_profiled = FALSE
      )
    }
  } else if (estmethod == "sv-wls") {
    # find sv-wls objective function value
    objval <- get_svloss(spcov_params_val,
      esv = esv,
      weights = weights
    )
  } else if (estmethod == "sv-cl") {
    # find cl objective function value
    objval <- get_glogclikloss(spcov_params_val, residual_vector2, dist_vector)
  }
  objval
}
