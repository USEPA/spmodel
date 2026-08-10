# shared building blocks for the cov_initial_search()/cov_initial_search_glm()
# grid-search family (R/cov_initial_search.R, R/cov_initial_search_glm.R).
# Extracted because the 4 covariance-type-group "body" functions in each file
# (exponential-family, none/ie, matern-family, car/sar), each branching into a
# no-randcov and a with-randcov path, duplicate several multi-line grid-
# construction blocks verbatim or near-verbatim -- see
# R/use_loglik_helpers.R/R/spcov_transform_helpers.R for the precedent this
# follows.

#' Build the initial spatial covariance parameter starting grid
#'
#' @param de Candidate de (spatially structured variance) proportions
#' @param ie Candidate ie (independent error variance) proportions
#' @param ns2 The inflated OLS sample variance (the total variance budget)
#' @param ... Additional named candidate-value vectors passed to
#'   \code{expand.grid()} (typically \code{range}, and \code{rotate}/
#'   \code{scale} under anisotropy or \code{extra} for matern-family types;
#'   car/sar pass only \code{range}, with no rotate/scale/extra columns)
#'
#' @return A data frame with one row per candidate combination where the
#'   de/ie proportions sum to 1 (i.e. together exhaust the variance budget
#'   with no double counting), with de/ie rescaled from proportions onto their
#'   final absolute-variance scale (\code{ns2 * proportion})
#'
#' @noRd
build_de_ie_grid <- function(de, ie, ns2, ...) {
  spcov_grid <- expand.grid(de = de, ie = ie, ...)
  spcov_grid <- spcov_grid[spcov_grid$de + spcov_grid$ie == 1, , drop = FALSE]
  if (all(c("rotate", "scale") %in% names(spcov_grid))) {
    # scale = 1 is a circle (no anisotropic stretching), and rotating a circle
    # leaves it unchanged, so every rotate candidate at scale = 1 evaluates an
    # identical covariance -- only rotate = 0 needs to be kept
    spcov_grid <- spcov_grid[spcov_grid$scale != 1 | spcov_grid$rotate == 0, , drop = FALSE]
  }
  spcov_grid[, c("de", "ie")] <- ns2 * spcov_grid[, c("de", "ie")]
  spcov_grid
}

#' Build and union the spatially-dominant / evenly-dominated / random-dominant
#' covariance parameter starting grids for models with random effects
#'
#' With random effects present, the total variance can plausibly be split many
#' ways between spatial dependence, nugget, and each random effect. Rather
#' than crossing every combination (which would blow up combinatorially),
#' three representative grids are built -- one where spatial variance
#' dominates, one where variance is spread evenly, and one where the random
#' effect(s) dominate -- and unioned together as the candidate starting
#' points.
#'
#' @param spcov_grid The initial spatial covariance starting grid from
#'   \code{build_de_ie_grid()}, on its final (already \code{ns2}-scaled) scale
#' @param spcov_grid_init A saved copy of \code{spcov_grid} from before any
#'   random-effect scaling was applied -- used to build the evenly-dominated
#'   and random-dominated grids from the original candidate combinations
#' @param ns2 The inflated OLS sample variance (the total variance budget)
#' @param spcov_initial_NA A spatial initial NA object
#' @param randcov_initial_NA A random effect initial NA object
#' @param nvar_spcov Number of spatial variance parameters (1 for none/ie,
#'   2 otherwise)
#' @param scale_cols Which \code{spcov_grid}/\code{spcov_grid_init} columns
#'   represent variance components to rescale at each stage
#'   (\code{c("de", "ie")}, or \code{c("de", "ie", "extra")} for car/sar,
#'   whose \code{extra} column is an alias of \code{de})
#' @param evencov_filter A logical vector (parallel to \code{spcov_grid_init}'s
#'   rows) selecting which rows feed the evenly-dominated grid (spread-out
#'   spatial cases); none/ie has no meaningful de/ie spread to filter on, so
#'   it defaults to keeping every row
#' @param randcov_filter Same, for the random-dominated grid's spatial
#'   component
#'
#' @return The unioned covariance parameter starting grid (spatially
#'   dominant + evenly dominated + random dominant), deduplicated
#'
#' @noRd
add_randcov_grids <- function(spcov_grid, spcov_grid_init, ns2, spcov_initial_NA, randcov_initial_NA,
                              nvar_spcov, scale_cols = c("de", "ie"),
                              evencov_filter = spcov_grid_init$de == spcov_grid_init$ie,
                              randcov_filter = spcov_grid_init$de == spcov_grid_init$ie &
                                spcov_grid_init$range == min(spcov_grid_init$range)) {
  randcov_names <- names(randcov_initial_NA$initial)
  nvar_randcov <- length(randcov_names)

  # spatially dominant grid
  ## spatial components 90% of variance
  spcov_grid[, scale_cols] <- 0.9 * spcov_grid[, scale_cols]

  # replace with initial
  for (x in names(spcov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }
  ## random effects 10% of variance and evenly spread
  for (x in randcov_names) {
    if (!is.na(randcov_initial_NA$initial[[x]])) {
      spcov_grid[, x] <- randcov_initial_NA$initial[[x]]
    } else {
      spcov_grid[, x] <- 0.1 * ns2 / nvar_randcov
    }
  }

  # Evenly dominated grid
  ## find all variance parameters
  nvar_cov <- nvar_spcov + nvar_randcov
  ## keep only spatial cases with variance spread out
  evencov_grid <- spcov_grid_init[evencov_filter, , drop = FALSE]
  ## scale to incorporate overall variance
  evencov_grid[, scale_cols] <- nvar_spcov / nvar_cov * evencov_grid[, scale_cols]

  # replace with initial values
  ## spatial
  for (x in names(evencov_grid)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      evencov_grid[, x] <- spcov_initial_NA$initial[[x]]
    }
  }
  ## random
  for (x in randcov_names) {
    if (!is.na(randcov_initial_NA$initial[[x]])) {
      evencov_grid[, x] <- randcov_initial_NA$initial[[x]]
    } else {
      evencov_grid[, x] <- 1 / nvar_cov * ns2
    }
  }
  # find unique combinations
  evencov_grid <- unique(evencov_grid)

  # random dominated grid
  ## spatial component 10%
  randcov_grid_spcov <- spcov_grid_init[randcov_filter, , drop = FALSE]
  randcov_grid_spcov[, scale_cols] <- 0.1 * randcov_grid_spcov[, scale_cols]

  # replace spatial initial values
  for (x in names(randcov_grid_spcov)) {
    if (!is.na(spcov_initial_NA$initial[[x]])) {
      randcov_grid_spcov[, x] <- spcov_initial_NA$initial[[x]]
    }
  }
  # find unique values
  randcov_grid_spcov <- unique(randcov_grid_spcov)

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
  for (x in randcov_names) {
    if (!is.na(randcov_initial_NA$initial[[x]])) {
      randcov_grid_randcov[, x] <- randcov_initial_NA$initial[[x]]
    }
  }

  ## bind together and replicate
  randcov_grid_spcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_randcov), randcov_grid_spcov, simplify = FALSE))
  randcov_grid_randcov_rep <- do.call("rbind", replicate(NROW(randcov_grid_spcov), randcov_grid_randcov, simplify = FALSE))
  randcov_grid <- cbind(randcov_grid_spcov_rep, randcov_grid_randcov_rep)

  # bind together all grids
  cov_grid <- rbind(spcov_grid, evencov_grid, randcov_grid)
  unique(cov_grid)
}

#' Add a dispersion candidate column to a GLM covariance parameter starting grid
#'
#' @param grid The covariance parameter starting grid
#' @param dispersion_initial_NA A dispersion initial NA object
#'
#' @return \code{grid} with a \code{dispersion} column added: 1 (i.e. no
#'   over/under-dispersion beyond what the GLM family implies) unless the
#'   user fixed a value, in which case that fixed value is used instead
#'
#' @noRd
add_dispersion_column <- function(grid, dispersion_initial_NA) {
  grid$dispersion <- 1
  if (!is.na(dispersion_initial_NA$initial)) {
    grid[, "dispersion"] <- dispersion_initial_NA$initial
  }
  grid
}

#' Precompute the empirical semivariogram used by the sv-wls objective
#'
#' The sv-wls objective (evaluated per grid point in \code{eval_grid()}) is a
#' weighted least-squares fit to an empirical semivariogram, so it only needs
#' to be computed once, up front, rather than inside the grid loop. A no-op
#' (returns \code{esv_val = NULL}) unless \code{estmethod} is \code{"sv-wls"}
#' -- Gaussian-only, since sv-wls isn't defined for GLM responses.
#'
#' @param estmethod The estimation method
#' @param data_object The data object
#' @param spcov_initial_NA A spatial initial NA object
#' @param dist_matrix_list A list of distance matrices
#' @param esv_dotlist Additional arguments passed to \code{esv()}
#'
#' @return A list with \code{esv_val} and \code{dist_matrix_list} (rotated
#'   under anisotropy to match \code{spcov_initial_NA}'s fixed rotate/scale,
#'   otherwise unchanged)
#'
#' @noRd
precompute_sv_wls <- function(estmethod, data_object, spcov_initial_NA, dist_matrix_list, esv_dotlist) {
  if (estmethod == "sv-wls") {
    if (data_object$anisotropy) {
      new_coords_list <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
        rotate = spcov_initial_NA$initial[["rotate"]],
        scale = spcov_initial_NA$initial[["scale"]]
      )
      dist_matrix_list <- lapply(new_coords_list, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))
    }
    # compute empirical semivariogram
    esv_vals <- mapply(d = data_object$obdata_list, m = dist_matrix_list, function(d, m) {
      do.call("esv", c(
        list(
          formula = data_object$formula,
          data = d,
          dist_matrix = m,
          partition_factor = data_object$partition_factor
        ),
        esv_dotlist
      ))
    }, SIMPLIFY = FALSE)
    esv_vals <- do.call("rbind", esv_vals)
    esv_vals <- esv_vals[esv_vals$np > 0, , drop = FALSE]
    esv_vals$bins <- droplevels(esv_vals$bins)
    esv_val <- data.frame(
      bins = levels(esv_vals$bins),
      dist = tapply(esv_vals$dist, esv_vals$bins, function(x) mean(x, na.rm = TRUE)),
      gamma = tapply(esv_vals$gamma, esv_vals$bins, function(x) mean(x, na.rm = TRUE)),
      np = tapply(esv_vals$np, esv_vals$bins, function(x) mean(x))
    )
  } else {
    esv_val <- NULL
  }
  list(esv_val = esv_val, dist_matrix_list = dist_matrix_list)
}

#' Precompute the pairwise distance/residual vectors used by the sv-cl objective
#'
#' For sv-cl (composite likelihood via pairwise squared differences, Curriero
#' & Lele 1999), the vector of pairwise distances and squared OLS-residual
#' differences is precomputed once, up front, since the composite-likelihood
#' objective in \code{eval_grid()} is a function of these fixed vectors, not
#' of the grid point itself. A no-op (returns \code{dist_vector =
#' residual_vector2 = NULL}) unless \code{estmethod} is \code{"sv-cl"} --
#' Gaussian-only, since sv-cl isn't defined for GLM responses.
#'
#' @param estmethod The estimation method
#' @param data_object The data object
#' @param spcov_initial_NA A spatial initial NA object
#' @param dist_matrix_list A list of distance matrices
#'
#' @return A list with \code{dist_vector}, \code{residual_vector2}, and
#'   \code{dist_matrix_list} (rotated under anisotropy to match
#'   \code{spcov_initial_NA}'s fixed rotate/scale, otherwise unchanged)
#'
#' @noRd
precompute_sv_cl <- function(estmethod, data_object, spcov_initial_NA, dist_matrix_list) {
  if (estmethod == "sv-cl") {
    if (data_object$anisotropy) {
      new_coords_list <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
        rotate = spcov_initial_NA$initial[["rotate"]],
        scale = spcov_initial_NA$initial[["scale"]]
      )
      dist_matrix_list <- lapply(new_coords_list, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))
    }
    dist_vector_list <- lapply(dist_matrix_list, function(x) {
      x <- as.matrix(x)
      # only the upper triangle is needed since distance/residual-difference
      # matrices are symmetric with a meaningless (zero) diagonal
      x <- x[upper.tri(x)]
    })
    residual_list <- lapply(data_object$obdata_list, function(d) residuals(lm(data_object$formula, data = d)))
    # spdist() on a single vector of residuals gives pairwise |residual_i - residual_j|
    residual_matrix_list <- lapply(residual_list, function(x) spdist(xcoord_val = x))
    residual_vector_list <- lapply(residual_matrix_list, function(x) {
      x <- as.matrix(x)
      x <- x[upper.tri(x)]
    })
    residual_vector <- unlist(residual_vector_list)
    if (!is.null(data_object$partition_list)) {
      # zero out (via multiplication by a 0/1 indicator) pairs that fall in
      # different partitions, since composite likelihood only uses within-partition pairs
      partition_vector_list <- lapply(data_object$partition_list, function(x) {
        x <- as.matrix(x)
        x <- x[upper.tri(x)]
      })
      dist_vector_list <- mapply(d = dist_vector_list, p = partition_vector_list, function(d, p) d * p, SIMPLIFY = FALSE)
      residual_vector_list <- mapply(r = residual_vector_list, p = partition_vector_list, function(r, p) r * p, SIMPLIFY = FALSE)
    }

    dist_vector <- unlist(dist_vector_list)
    # drop zero distances (self-pairs, or pairs zeroed out by partitioning above)
    dist_index <- dist_vector > 0
    dist_vector <- dist_vector[dist_index]
    residual_vector <- unlist(residual_vector_list)
    residual_vector <- residual_vector[dist_index]
    residual_vector2 <- residual_vector^2
  } else {
    dist_vector <- NULL
    residual_vector2 <- NULL
  }
  list(dist_vector = dist_vector, residual_vector2 = residual_vector2, dist_matrix_list = dist_matrix_list)
}

#' Evaluate every grid point and select the winning combination
#'
#' @param grid The covariance parameter starting grid
#' @param eval_fn The per-grid-point objective function (\code{eval_grid()}
#'   or \code{eval_grid_glm()})
#' @param ... Additional named arguments forwarded to \code{eval_fn} via
#'   \code{vapply()} -- these differ by family: \code{eval_grid()} takes
#'   \code{weights}/\code{esv}/\code{dist_vector}/\code{residual_vector2};
#'   \code{eval_grid_glm()} takes \code{family} instead
#'
#' @return The winning grid row, unlisted into a named numeric vector (as
#'   callers expect, e.g. \code{min_params[["de"]]})
#'
#' @noRd
select_best_grid_point <- function(grid, eval_fn, ...) {
  grid_splits <- split(grid, seq_len(NROW(grid)))
  objvals <- vapply(X = grid_splits, FUN = eval_fn, FUN.VALUE = numeric(1), ...)
  unlist(grid_splits[[which.min(objvals)]])
}

#' Resolve the anisotropy rotation ambiguity when evaluating one grid point
#'
#' The rotate angle is only identifiable modulo pi (rotating by \code{rotate}
#' vs. \code{pi - rotate} gives the same ellipse orientation), so the
#' objective is evaluated at both candidate angles and the smaller (better)
#' value is kept. Scoped locally to \code{eval_grid()}/\code{eval_grid_glm()}
#' -- unlike \code{resolve_anis_rotation()}, which also needs to
#' report which candidate wins and its \code{dist_matrix_list} once
#' \code{optim()} has converged, a single grid-point evaluation only ever
#' needs the scalar minimum, so this is kept as a separate, simpler helper
#' rather than shared across both call sites.
#'
#' @param spcov_params_val A \code{spcov_params} object
#' @param randcov_params_val A \code{randcov_params} object, or \code{NULL}
#' @param data_object The data object
#' @param estmethod The estimation method
#' @param dispersion_params_val A \code{dispersion_params} object, or
#'   \code{NULL} for Gaussian (non-GLM) estimation -- when supplied,
#'   dispatches to the Laplace-approximation products/loss functions instead
#'   of the Gaussian ones
#'
#' @return The smaller (better) of the two candidate angles' \code{-2} times
#'   log-likelihood values
#'
#' @noRd
resolve_rotation_ambiguity <- function(spcov_params_val, randcov_params_val, data_object, estmethod,
                                       dispersion_params_val = NULL) {
  new_coords_list_q1 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
    rotate = spcov_params_val[["rotate"]], scale = spcov_params_val[["scale"]]
  )
  dist_matrix_list_q1 <- lapply(new_coords_list_q1, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

  new_coords_list_q2 <- lapply(data_object$obdata_list, transform_anis, data_object$xcoord, data_object$ycoord,
    rotate = abs(pi - spcov_params_val[["rotate"]]), scale = spcov_params_val[["scale"]]
  )
  dist_matrix_list_q2 <- lapply(new_coords_list_q2, function(x) spdist(xcoord_val = x$xcoord_val, ycoord_val = x$ycoord_val))

  if (is.null(dispersion_params_val)) {
    prods_q1 <- gloglik_products(spcov_params_val, data_object, estmethod, dist_matrix_list_q1, randcov_params_val)
    prods_q2 <- gloglik_products(spcov_params_val, data_object, estmethod, dist_matrix_list_q2, randcov_params_val)
    objval_q1 <- get_minustwologlik(prods_q1, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE, randcov_profiled = FALSE)
    objval_q2 <- get_minustwologlik(prods_q2, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE, randcov_profiled = FALSE)
  } else {
    prods_q1 <- laploglik_products(spcov_params_val, dispersion_params_val, data_object, estmethod, dist_matrix_list_q1, randcov_params_val)
    prods_q2 <- laploglik_products(spcov_params_val, dispersion_params_val, data_object, estmethod, dist_matrix_list_q2, randcov_params_val)
    objval_q1 <- get_minustwolaploglik(prods_q1, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE, randcov_profiled = FALSE)
    objval_q2 <- get_minustwolaploglik(prods_q2, estmethod, data_object$n, data_object$p, spcov_profiled = FALSE, randcov_profiled = FALSE)
  }

  min(objval_q1, objval_q2)
}
