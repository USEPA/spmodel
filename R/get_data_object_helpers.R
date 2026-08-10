#' Reject legacy sp objects in favor of sf
#'
#' @param data data
#'
#' @return Error message or nothing
#'
#' @noRd
check_sp_not_supported <- function(data) {
  # the legacy sp package is not supported directly -- detect it by checking
  # which package the data's class attribute came from and error out early
  attr_sp <- attr(class(data), "package")
  if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
    stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
  }
  invisible(NULL)
}

#' Detect sf metadata for areal data, without any point-referenced coordinate handling
#'
#' @param data data
#'
#' @return A list with elements \code{is_sf}, \code{sf_column_name}, and \code{crs}
#'
#' @noRd
get_areal_sf_info <- function(data) {
  if (inherits(data, "sf")) {
    is_sf <- TRUE
    sf_column_name <- attributes(data)$sf_column
    crs <- attributes(data[[sf_column_name]])$crs
  } else {
    is_sf <- FALSE
    sf_column_name <- NULL
    crs <- NULL
  }
  list(is_sf = is_sf, sf_column_name = sf_column_name, crs = crs)
}

#' Expand a "." in formula into explicit predictor names, excluding
#' reserved coordinate/geometry columns
#'
#' @param formula A two-sided formula, possibly containing \code{.} on the
#'   right-hand side
#' @param data The data \code{.} is expanded against (only its column names
#'   matter here, not its values, so any data frame with the right columns
#'   works -- e.g. \code{obdata} or the full \code{data})
#' @param reserved_cols Column names to exclude from \code{.}'s expansion:
#'   \code{c(xcoord, ycoord, ycoord_orig_name)} for \code{splm()}/\code{spglm()},
#'   or the sf geometry column name for \code{spautor()}/\code{spgautor()}
#'   (\code{character(0)} if there is nothing to reserve, e.g. non-sf areal data)
#'
#' @details If \code{formula} actually contains \code{.},
#'   \code{.} is resolved via \code{terms()}, given a zero-row plain data frame
#'   containing only \code{data}'s non-reserved column \emph{names} (not a
#'   subset of \code{data} itself -- some data classes' \code{[} methods don't
#'   honor column exclusion the way a plain data frame's does; notably, an
#'   \code{sf} object's \code{[} always keeps the geometry column even when
#'   it isn't selected). Since \code{terms()} only consults a \code{data}
#'   argument's column names to decide what \code{.} stands for (not to
#'   validate other, explicitly-named terms), this leaves any explicit use of
#'   a reserved column elsewhere in \code{formula} (e.g. a trend-surface term
#'   deliberately using \code{xcoord} as a covariate) untouched, and the
#'   response is automatically excluded from \code{.} by R's usual formula
#'   semantics. The returned formula must be used in place of \code{formula}
#'   for the rest of model fitting/prediction (not just the one
#'   \code{model.frame()} call it's needed for) -- \code{object$formula} is
#'   reused verbatim by \code{predict()}, \code{kcv()}/\code{loocv()} refits,
#'   \code{anova()}, and \code{model.matrix()}, and if \code{.} were left
#'   unexpanded there, each of those would silently re-expand it against
#'   whatever columns their own data happens to have (e.g. \code{newdata}),
#'   which could pull the coordinate columns back in as predictors or
#'   otherwise produce a different predictor set than the one actually fit.
#'
#' @return \code{formula}, expanded if it contained \code{.}
#'
#' @noRd
expand_formula_dot <- function(formula, data, reserved_cols = character(0)) {
  if (!"." %in% all.vars(formula)) {
    return(formula)
  }
  dot_names <- setdiff(names(data), reserved_cols)
  dot_data <- as.data.frame(matrix(nrow = 0, ncol = length(dot_names), dimnames = list(NULL, dot_names)))
  formula(terms(formula, data = dot_data))
}

#' Resolve point-referenced coordinates: sf centroid coercion, xcoord/ycoord
#' validation, and \code{dim_coords} derivation
#'
#' @param data data
#' @param spcov_initial A \code{spcov_initial} object
#' @param estmethod The estimation method
#' @param xcoord The x-coordinate name, already resolved by the caller to
#'   either a plain string or \code{NULL} (if the original argument was
#'   omitted) -- see \code{@details}
#' @param ycoord The y-coordinate name, already resolved the same way
#'
#' @details Callers must resolve \code{xcoord}/\code{ycoord} before calling
#'   this function: \code{if (missing(xcoord)) NULL else
#'   as.character(substitute(xcoord))}, checked at the point where the
#'   original, possibly-omitted argument was captured (the one place
#'   \code{missing()} can answer that reliably). \code{NULL} signals
#'   "omitted" -- checked here via \code{is.null(xcoord)}, not
#'   \code{missing()}, since \code{xcoord}/\code{ycoord} are always actually
#'   supplied to this function (never R-missing) by the time they get here.
#'
#' @return A list with elements \code{data} (possibly converted from sf to a
#'   data frame, and/or with dummy coordinate columns added), \code{xcoord},
#'   \code{ycoord}, \code{dim_coords}, \code{ycoord_orig_name},
#'   \code{ycoord_orig_val}, \code{is_sf}, \code{sf_column_name}, \code{crs},
#'   and \code{data_sf}
#'
#' @noRd
get_point_ref_coords <- function(data, spcov_initial, estmethod, xcoord, ycoord) {
  # convert sf to data frame (point geometry) (1d objects obsolete)
  ## see if data has sf class
  if (inherits(data, "sf")) {
    # set is_sf
    is_sf <- TRUE
    sf_column_name <- attributes(data)$sf_column
    crs <- attributes(data[[sf_column_name]])$crs
    if (!inherits(spcov_initial, c("none", "ie")) && any(sf::st_geometry_type(data) != "POINT")) {
      warning("At least one geometry type in data is not equal to \"POINT\". Attempting to coerce all non-\"POINT\" geometries to \"POINT\" geometries via their centroids using sf::st_centroid().", call. = FALSE)
    }
    # all downstream distance/covariance computations need point coordinates,
    # so collapse any polygon/line geometries to their centroids
    data_sf <- suppressWarnings(sf::st_centroid(data))
    # store as data frame
    data <- sf_to_df(data_sf)
    if (!is.null(xcoord) || !is.null(ycoord)) {
      warning("data is an sf object. Ignoring xcoord and ycoord arguments.", call. = FALSE)
    }
    ## name xcoord ".xcoord" to be used later
    xcoord <- ".xcoord"
    ## name ycoord ".ycoord" to be used later
    ycoord <- ".ycoord"
  } else {
    is_sf <- FALSE
    sf_column_name <- NULL
    crs <- NULL
    data_sf <- NULL
  }

  if (!is_sf && is.null(xcoord) && !inherits(spcov_initial, c("none", "ie"))) {
    stop("The xcoord argument must be specified.", call. = FALSE)
  }

  if (!is.null(xcoord)) {
    if (!xcoord %in% colnames(data)) {
      stop("The xcoord argument must match the name of a variable in data.", call. = FALSE)
    }
  }

  if (!is.null(ycoord)) {
    if (!ycoord %in% colnames(data)) {
      stop("The ycoord argument must match the name of a variable in data.", call. = FALSE)
    }
  }

  # setting ycoord orig val for use with circular or triangular
  ycoord_orig_name <- NULL
  ycoord_orig_val <- NULL
  # find coordinate dimension and set defaults
  # dim_coords tracks how many spatial dimensions actually matter for this
  # model: 0 when there is no spatial covariance structure at all (so
  # coordinates are irrelevant and dummied out to zero), 1 when only a single
  # coordinate is used (e.g. triangular/cosine covariance, or ycoord omitted),
  # and 2 for standard two-dimensional spatial covariances
  if (inherits(spcov_initial, c("none", "ie")) && estmethod %in% c("reml", "ml")) {
    dim_coords <- 0
    if (is.null(xcoord)) {
      xcoord <- ".xcoord"
      data[[xcoord]] <- 0
    }
    if (is.null(ycoord)) {
      ycoord <- ".ycoord"
      if (xcoord == ".ycoord") {
        ycoord <- ".ycoord2"
      }
      data[[ycoord]] <- 0
    }
  } else if (is.null(ycoord) || inherits(spcov_initial, c("triangular", "cosine"))) {
    dim_coords <- 1
    # if the user did supply ycoord but the covariance type only uses one
    # dimension (triangular/cosine), keep the original values around so they
    # can still be referenced/restored later even though they're not used here
    if (!is.null(ycoord)) {
      ycoord_orig_name <- ycoord
      ycoord_orig_val <- data[[ycoord]]
    }
    ycoord <- ".ycoord"
    if (xcoord == ".ycoord") {
      ycoord <- ".ycoord2"
    }
    data[[ycoord]] <- 0
  } else {
    dim_coords <- 2
  }

  # check missing coordinates (missing coordinates can't be in sf objects)
  if (any(is.na(c(data[[xcoord]], data[[ycoord]])))) {
    stop("Missing values in coordinates not allowed.", call. = FALSE)
  }

  # check coordinates proper type
  if (any(!is.numeric(data[[xcoord]]), !is.numeric(data[[ycoord]]))) {
    stop("Coordinates must be numeric.", call. = FALSE)
  }

  # check if coordinates are projected
  if (is_sf) {
    if (!is.na(st_is_longlat(crs)) && st_is_longlat(crs)) {
      warning("Coordinates are in a geographic coordinate system and will be used as is. For the most accurate results, please ensure coordinates are in a projected coordinate system (e.g., via sf::st_transform()).", call. = FALSE)
    }
  }

  list(
    data = data, xcoord = xcoord, ycoord = ycoord, dim_coords = dim_coords,
    ycoord_orig_name = ycoord_orig_name, ycoord_orig_val = ycoord_orig_val,
    is_sf = is_sf, sf_column_name = sf_column_name, crs = crs, data_sf = data_sf
  )
}

#' Compute the range-parameter constraint bound from the observed coordinates'
#' bounding-box half-diagonal
#'
#' @param obdata The observed data
#' @param xcoord The x-coordinate name
#' @param ycoord The y-coordinate name
#' @param spcov_initial A \code{spcov_initial} object
#' @param range_constrain Whether to constrain the range parameter
#'
#' @return A list with elements \code{max_halfdist}, \code{range_constrain},
#'   and \code{range_constrain_value}
#'
#' @noRd
get_range_constrain_setup <- function(obdata, xcoord, ycoord, spcov_initial, range_constrain) {
  # half the diagonal (largest possible pairwise distance) of the bounding
  # box of observed coordinates -- used below to put a sensible upper bound
  # on the range parameter, since ranges much larger than the domain itself
  # are not identifiable from the data
  x_range <- range(obdata[[xcoord]])
  y_range <- range(obdata[[ycoord]])
  max_halfdist <- sqrt((max(x_range) - min(x_range))^2 + (max(y_range) - min(y_range))^2) / 2

  max_range_scale <- 4
  range_constrain_value <- 2 * max_halfdist * max_range_scale
  # skip constraining if the range is already fixed (known) or if the user's
  # own initial range guess already exceeds the constraint bound
  if ("range" %in% names(spcov_initial$is_known)) {
    if (spcov_initial$is_known[["range"]] || (spcov_initial$initial[["range"]] > range_constrain_value)) {
      range_constrain <- FALSE
    }
  }

  if (inherits(spcov_initial, c("none", "ie"))) {
    range_constrain <- FALSE
  }

  if (is.logical(range_constrain)) {
    if (!range_constrain) {
      range_constrain_value <- NULL
    }
  } else {
    stop("range_constrain must be logical.", call. = FALSE)
  }

  list(max_halfdist = max_halfdist, range_constrain = range_constrain, range_constrain_value = range_constrain_value)
}

#' Coerce a partition factor formula to its canonical reformulated form
#'
#' @param partition_factor A partition factor formula (or \code{NULL})
#' @param obdata The observed data
#'
#' @return The reformulated \code{partition_factor}, or \code{NULL}
#'
#' @noRd
coerce_partition_factor <- function(partition_factor, obdata) {
  if (!is.null(partition_factor)) {
    partition_factor_labels <- labels(terms(partition_factor))
    if (length(partition_factor_labels) > 1) {
      stop("Only one variable can be specified in partition_factor.", call. = FALSE)
    }
    partition_mf <- model.frame(partition_factor, obdata)
    if (any(!attr(terms(partition_mf), "dataClasses") %in% c("character", "factor", "ordered"))) {
      stop("Partition factor variable must be categorical or factor.", call. = FALSE)
    }
    partition_factor <- reformulate(partition_factor_labels, intercept = FALSE)
  }
  partition_factor
}

#' Validate and rename a \code{randcov_initial} object's parameter names
#'
#' @param randcov_initial A \code{randcov_initial} object (or \code{NULL})
#' @param dedupe_names Whether to \code{unique()} the derived names before
#'   comparing lengths (point-referenced constructors do; areal constructors
#'   do not -- a pre-existing difference, preserved here rather than
#'   harmonized)
#'
#' @return A validated (and possibly freshly-created) \code{randcov_initial} object
#'
#' @noRd
validate_randcov_initial <- function(randcov_initial, dedupe_names = TRUE) {
  if (is.null(randcov_initial)) {
    randcov_initial <- spmodel::randcov_initial()
  } else {
    randcov_given_names <- unlist(lapply(
      names(randcov_initial$initial),
      function(x) labels(terms(reformulate(x)))
    ))
    randcov_initial_names <- if (dedupe_names) {
      unique(unlist(lapply(randcov_given_names, get_randcov_name)))
    } else {
      unlist(lapply(randcov_given_names, get_randcov_name))
    }
    if (length(randcov_initial_names) != length(names(randcov_initial$initial))) {
      stop("No / can be specified in randcov_initial(). Please specify starting
             values for each variable (e.g., a/b = a + a:b)", call. = FALSE)
    }
    names(randcov_initial$initial) <- randcov_initial_names
    names(randcov_initial$is_known) <- randcov_initial_names
  }
  randcov_initial
}

#' Require the response to be numeric with nonzero variance
#'
#' @param y The response vector
#'
#' @return Error message or nothing
#'
#' @noRd
check_response_numeric_and_variable <- function(y) {
  if (!is.numeric(y)) {
    stop("Response variable must be numeric", call. = FALSE)
  }
  if (var(y) == 0) {
    stop("The response has no variability. Model fit unreliable.", call. = FALSE)
  }
  invisible(NULL)
}

#' Warn about perfect collinearity in the design matrix
#'
#' @param p The rank of \code{X}
#' @param X The design matrix
#'
#' @return Warning or nothing
#'
#' @noRd
check_rank_collinearity <- function(p, X) {
  if (p < NCOL(X)) {
    warning("There are perfect collinearities detected in X (the matrix of explanatory variables). This may make the model fit unreliable or may cause an error while model fitting. Consider removing redundant explanatory variables and refitting the model.", call. = FALSE)
  }
  invisible(NULL)
}

#' Require the number of fixed effects to be less than the sample size
#'
#' @param p The rank of \code{X}
#' @param n The sample size
#' @param refit_fun_text The function name to suggest rerunning, used
#'   literally in the error message (some existing call sites reference a
#'   different constructor's name than their own -- a pre-existing
#'   inconsistency, preserved here rather than harmonized)
#'
#' @return Error message or nothing
#'
#' @noRd
check_p_n <- function(p, n, refit_fun_text) {
  if (p >= n) {
    stop(
      paste0(
        "The number of fixed effects is at least as large as the number of observations (p >= n). Consider reducing the number of fixed effects and rerunning ",
        refit_fun_text, "()."
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Require a user-supplied local$index to match the length of the non-missing response
#'
#' @param local The (possibly list) \code{local} argument
#' @param n_obdata The number of rows in \code{obdata} (the non-missing-response,
#'   complete-predictor data actually used for fitting)
#'
#' @details \code{local$index} labels each row of the data used for fitting
#'   with a partition assignment (see the spatial indexing/SPIN method this
#'   implements), so its length must match \code{obdata}, not the original
#'   (possibly larger) \code{data} passed to \code{splm()}/\code{spglm()} --
#'   rows with a missing response are excluded from fitting entirely (they
#'   become prediction locations instead) and so must already be excluded
#'   from \code{local$index} too
#'
#' @return Error message or nothing
#'
#' @noRd
check_local_index_length <- function(local, n_obdata) {
  if (is.list(local) && "index" %in% names(local)) {
    if (length(local$index) != n_obdata) {
      stop(
        "local$index must have the same length as the non-missing (non-NA) response vector (",
        n_obdata, "), but has length ", length(local$index), ". ",
        "Observations with a missing response are excluded from fitting (they become prediction locations instead), so they must also be excluded from local$index.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

#' Resolve the local (big-data) partitioning of the observed data, once
#' \code{X}/\code{y}/\code{offset}/\code{size} are known
#'
#' @param local A list of big-data options (or \code{NULL})
#' @param obdata The observed data
#' @param xcoord The x-coordinate name
#' @param ycoord The y-coordinate name
#' @param n The sample size
#' @param partition_factor A partition factor formula (or \code{NULL})
#' @param n_threshold The sample-size threshold above which \code{local} auto-triggers
#' @param local_message The exact message to emit when auto-triggering (callers'
#'   wording differs by their own threshold value)
#' @param X The design matrix
#' @param y The response vector
#' @param offset The offset (or \code{NULL})
#' @param size Binomial trial counts (or \code{NULL}, GLM constructors only)
#'
#' @return A list with elements \code{local}, \code{obdata_list}, \code{X_list},
#'   \code{y_list}, \code{ones_list}, \code{offset}, and \code{size}
#'
#' @noRd
build_local_partition_lists <- function(local, obdata, xcoord, ycoord, n, partition_factor,
                                         n_threshold, local_message, X, y, offset = NULL, size = NULL) {
  # exact likelihood evaluation requires an n x n covariance matrix, which is
  # infeasible for large n; above this threshold, default to partitioning the
  # data into local neighborhoods for a computationally cheaper approximation
  if (is.null(local)) {
    if (n > n_threshold) {
      local <- TRUE
      message(local_message)
    } else {
      local <- FALSE
    }
  }
  local <- get_local_list_estimation(local, obdata, xcoord, ycoord, n, partition_factor)

  # store data list
  obdata_list <- split.data.frame(obdata, local$index)

  # store X and y
  X_list <- split.data.frame(X, local$index)
  y_list <- split.data.frame(y, local$index)
  ones_list <- lapply(obdata_list, function(x) matrix(rep(1, nrow(x)), ncol = 1))

  if (!is.null(size)) {
    # size must be reordered/split the same way as y so trial counts stay
    # aligned with their observations after partitioning
    size_list <- split(size, local$index) # just split because vector not matrix
    size <- as.vector(do.call("c", size_list)) # rearranging size by y list
  }

  # organize offset (as a one col matrix)
  if (!is.null(offset)) {
    offset <- do.call("rbind", (split.data.frame(offset, local$index)))
  }

  list(
    local = local, obdata_list = obdata_list, X_list = X_list, y_list = y_list,
    ones_list = ones_list, offset = offset, size = size
  )
}

#' Build the CAR/SAR neighbor structure (\code{W}, \code{M}, and their derived
#' quantities) shared by \code{spautor()} and \code{spgautor()}
#'
#' @param data data
#' @param spcov_initial A \code{spcov_initial} object
#' @param W A neighbor weight matrix (or \code{NULL} to build one internally)
#' @param M A diagonal weighting matrix for the CAR covariance (or \code{NULL})
#' @param row_st Whether to row-standardize \code{W}
#' @param range_positive Whether the range parameter is constrained positive
#' @param cutoff The neighbor cutoff distance (when \code{W} is built internally)
#'
#' @return A list with elements \code{W}, \code{W_rowsums}, \code{is_W_connected},
#'   \code{M}, \code{rho_lb}, and \code{rho_ub}
#'
#' @noRd
build_car_neighbor_structure <- function(data, spcov_initial, W, M, row_st, range_positive, cutoff) {
  # create distance matrix (if not provided) -- sf::st_intersects() assumes
  # units are neighbors with themselves, so we need to set the diagonal of the
  # matrix equal to zero
  # W is the binary neighbor (adjacency) matrix for the areal/autoregressive
  # model: for point geometries, two sites are neighbors if they are within
  # cutoff distance of each other; for polygons, neighbors are sites whose
  # boundaries physically touch/overlap (st_intersects)
  if (is.null(W)) {
    geom_type <- st_geometry_type(data, by_geometry = FALSE)
    if (geom_type == "POINT") {
      if (is.null(cutoff)) {
        stop("cutoff must be specified if using a distance-based neighbor cutoff.", call. = FALSE)
      }
      coords_val <- st_coordinates(data)
      W <- 1 * (as.matrix(dist(coords_val)) <= cutoff)
      diag(W) <- 0
      if (sum(W) == 0) {
        stop("cutoff must be larger than the smallest distance between potential neighbors.", call. = FALSE)
      }
    } else {
      W <- sf::st_intersects(data, sparse = FALSE)
      diag(W) <- 0
    }
  }

  # turn W into a sparse Matrix and logical regardless of whether provided by us or user
  W <- 1 * Matrix::Matrix(W, sparse = TRUE)
  W_rowsums <- Matrix::rowSums(W)
  # sites with a zero row sum have no neighbors; CAR/SAR models generally
  # assume every site is connected to at least one other, so this flag is
  # checked elsewhere to warn/error appropriately
  is_W_connected <- all(W_rowsums > 0)

  # make M if necessary
  # M holds the (diagonal) conditional-variance weights for the CAR model;
  # row-standardizing W implies M = 1/rowSums(W) so that (I - rho*W) stays
  # symmetric relative to M, which the CAR formulation requires
  if (row_st) {
    if (!is.null(M)) {
      if (inherits(spcov_initial, "car")) {
        warning("Overriding M when row_st = TRUE", call. = FALSE)
      }
      if (inherits(spcov_initial, "sar")) {
        warning("M ignored for sar models", call. = FALSE)
      }
    }
    M <- 1 / W_rowsums # this has not been standardized
  } else {
    if (is.null(M)) {
      M <- rep(1, nrow(W)) # assume identity
    } else {
      if (inherits(spcov_initial, "sar")) {
        warning("M ignored for sar models", call. = FALSE)
      }
      M <- as.matrix(M) # coerce to matrix from vector, matrix, or Matrix
      if (dim(M)[1] == dim(M)[2]) {
        M <- diag(M) # take diagonal of matrix
      } else {
        M <- as.vector(M) # assume diagonal already given as one-column vector
      }
    }
  }

  # row standardize W if necessary
  if (row_st) {
    W_rowsums_val <- W_rowsums # make copy so rowsums are saved later
    W_rowsums_val[W_rowsums_val == 0] <- 1 # not a Matrix object so this subsetting is okay
    W <- W / W_rowsums_val
  }

  # the CAR covariance is proportional to (I - rho*W)^{-1} %*% diag(M); for
  # this to be a valid (symmetric, positive-definite-able) covariance matrix,
  # M^{-1}(I - rho*W) must be symmetric, which is checked here
  if (inherits(spcov_initial, "car") && !isSymmetric(as.matrix((Matrix(diag(nrow(W)), sparse = TRUE) - W) * 1 / M))) {
    stop("W and M must satisfy the CAR symmetry condition", call. = FALSE)
  }

  # find eigenvalues of W for connected sites
  # (I - rho*W) is invertible only for rho strictly between 1/min(eigenvalue)
  # and 1/max(eigenvalue), so these eigenvalues bound the feasible range for
  # the spatial autocorrelation parameter rho used during optimization
  rowsums_nonzero <- which(W_rowsums != 0)
  W_eigen <- Re(eigen(W[rowsums_nonzero, rowsums_nonzero])$values)
  if (range_positive) {
    rho_lb <- 1e-5
  } else {
    rho_lb <- 1 / min(W_eigen) + 1e-5 # rho strictly > lb
  }
  rho_ub <- 1 / max(W_eigen) - 1e-5 # rho strictly < ub

  list(W = W, W_rowsums = W_rowsums, is_W_connected = is_W_connected, M = M, rho_lb = rho_lb, rho_ub = rho_ub)
}
