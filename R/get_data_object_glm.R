#' Build the internal data object used to fit an \code{spglm()} model
#'
#' @param formula A formula
#' @param family The response family
#' @param data A data frame or \code{sf} object
#' @param spcov_initial A \code{spcov_initial} object
#' @param xcoord The x-coordinate name
#' @param ycoord The y-coordinate name
#' @param estmethod The estimation method
#' @param anisotropy Whether to model anisotropy
#' @param random A random effect formula (or \code{NULL})
#' @param randcov_initial A \code{randcov_initial} object (or \code{NULL})
#' @param partition_factor A partition factor formula (or \code{NULL})
#' @param local A list of big-data options (or \code{NULL})
#' @param range_constrain Whether to constrain the range parameter
#' @param ... Additional arguments
#'
#' @return A data object containing the (possibly partitioned) design matrices,
#'   response vector, coordinates, and other quantities needed throughout
#'   model fitting and prediction
#'
#' @noRd
get_data_object_spglm <- function(formula, family, data, spcov_initial, xcoord, ycoord, estmethod,
                                  anisotropy, random, randcov_initial, partition_factor, local,
                                  range_constrain, ...) {
  check_sp_not_supported(data)

  point_ref <- get_point_ref_coords(data, spcov_initial, estmethod, xcoord = xcoord, ycoord = ycoord)
  data <- point_ref$data
  xcoord <- point_ref$xcoord
  ycoord <- point_ref$ycoord
  dim_coords <- point_ref$dim_coords
  ycoord_orig_name <- point_ref$ycoord_orig_name
  ycoord_orig_val <- point_ref$ycoord_orig_val
  is_sf <- point_ref$is_sf
  sf_column_name <- point_ref$sf_column_name
  crs <- point_ref$crs
  data_sf <- point_ref$data_sf

  # expanding "." in formula
  formula <- expand_formula_dot(formula, data, c(xcoord, ycoord, ycoord_orig_name))

  # subsetting by na and not na values
  ## find response variable name
  # a binomial response can be specified as cbind(successes, failures), which
  # model.response() returns as a 2-column matrix; rowSums() collapses that
  # to a single NA-detection vector so both response encodings are handled
  # rows with a missing response are later folded into newdata for prediction
  response_index <- model.response(model.frame(formula, data, na.action = na.pass))
  if (NCOL(response_index) == 2) {
    na_index <- is.na(rowSums(response_index)) # will be NA if successes or failures NA
  } else {
    na_index <- is.na(response_index)
  }

  # store observed index
  observed_index <- which(!na_index)
  missing_index <- which(na_index)

  # find small and newdata
  if (any(na_index)) {
    ## find newdata to be used in prediction later
    if (is_sf) {
      newdata <- data_sf[na_index, , drop = FALSE] # keep as sf object is users want that
    } else {
      newdata <- data[na_index, , drop = FALSE]
    }
    ## subset original data
    obdata <- data[!na_index, , drop = FALSE]
  } else {
    obdata <- data
    newdata <- NULL
  }

  # finding model frame
  obdata_model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass)
  # finding contrasts as ...
  dots <- list(...)
  if (!"contrasts" %in% names(dots)) {
    dots$contrasts <- NULL
  }
  # model matrix with potential NA
  X <- model.matrix(formula, obdata_model_frame, contrasts = dots$contrasts)
  # finding rows w/out NA
  # unlike a missing response, missing predictors cannot become prediction
  # sites (predictions still require complete predictors), so this errors
  ob_predictors <- complete.cases(X)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }
  # subset obdata by nonNA predictors
  obdata <- obdata[ob_predictors, , drop = FALSE]
  # a user-supplied local$index must already be sized to obdata (see
  # check_local_index_length())
  check_local_index_length(local, NROW(obdata))

  # new model frame
  obdata_model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.omit)
  # find terms
  terms_val <- terms(obdata_model_frame)
  # find X
  X <- model.matrix(formula, obdata_model_frame, contrasts = dots$contrasts)
  # find induced contrasts and xlevels
  dots$contrasts <- attr(X, "contrasts")
  xlevels <- .getXlevels(terms_val, obdata_model_frame)
  # find p
  p <- as.numeric(Matrix::rankMatrix(X, method = "qr"))
  check_rank_collinearity(p, X)
  # find sample size
  n <- NROW(X)
  # find response
  # size is the number of binomial trials per observation, needed to convert
  # counts to proportions during model fitting; it is derived differently
  # depending on how the binomial response was specified by the user
  y_modr <- model.response(obdata_model_frame)
  if (NCOL(y_modr) == 2) {
    # cbind(successes, failures) form: y is successes, size is total trials
    y <- y_modr[, 1, drop = FALSE]
    size <- rowSums(y_modr)
  } else {
    if (family == "binomial") {
      # a two-level factor or logical response is a single-trial (0/1)
      # binomial outcome; coerce to numeric 0/1 so downstream math works
      if (is.factor(y_modr)) {
        if (length(levels(y_modr)) != 2) {
          stop("When family is binomial, a factor response must have exactly two levels.", call. = FALSE)
        }
        y_modr <- ifelse(y_modr == levels(y_modr)[1], 0, 1)
      }
      if (is.logical(y_modr)) {
        y_modr <- ifelse(y_modr, 1, 0) # or as.numeric()
      }
      size <- rep(1, n)
    } else {
      size <- NULL
    }
    y <- as.matrix(y_modr, ncol = 1)
  }

  # handle offset
  offset <- model.offset(obdata_model_frame)
  if (!is.null(offset)) {
    offset <- as.matrix(offset, ncol = 1)
  }

  check_response_numeric_and_variable(y)

  # checks on y
  response_checks_glm(family, y, size)

  check_p_n(p, n, "splm")

  # find s2 for initial values
  # log(y + 1) roughly approximates the link-scale response for common GLM
  # families (log link); an OLS fit on this transformed scale gives a rough
  # residual variance used only as an optimizer starting value, not the fit
  y_trans <- log(y + 1)
  qr_val <- qr(X)
  R_val <- qr.R(qr_val)
  betahat <- backsolve(R_val, qr.qty(qr_val, y_trans))
  resid <- y_trans - X %*% betahat
  s2 <- sum(resid^2) / (n - p)
  # small positive nugget added to the covariance diagonal for numerical
  # stability during GLM likelihood optimization (Laplace/other approximations)
  diagtol <- 1e-4

  range_setup <- get_range_constrain_setup(obdata, xcoord, ycoord, spcov_initial, range_constrain)
  max_halfdist <- range_setup$max_halfdist
  range_constrain <- range_setup$range_constrain
  range_constrain_value <- range_setup$range_constrain_value

  # override anisotropy argument if needed
  anisotropy <- get_anisotropy_corrected(anisotropy, spcov_initial)

  # coerce to factor
  partition_factor <- coerce_partition_factor(partition_factor, obdata)

  local_setup <- build_local_partition_lists(
    local, obdata, xcoord, ycoord, n, partition_factor,
    n_threshold = 3000,
    local_message = "Because the sample size exceeds 3,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.",
    X = X, y = y, offset = offset, size = size
  )
  local <- local_setup$local
  obdata_list <- local_setup$obdata_list
  X_list <- local_setup$X_list
  y_list <- local_setup$y_list
  ones_list <- local_setup$ones_list
  offset <- local_setup$offset
  size <- local_setup$size

  # store random effects list
  if (is.null(random)) {
    randcov_initial <- NULL
    randcov_list <- NULL
    randcov_names <- NULL
  } else {
    randcov_names <- get_randcov_names(random)
    randcov_xlevs <- lapply(randcov_names, get_randcov_xlev, obdata)
    names(randcov_xlevs) <- randcov_names
    randcov_list <- lapply(obdata_list, function(x) {
      get_randcov_Zs(x, randcov_names, xlev_list = randcov_xlevs)
    })
    randcov_initial <- validate_randcov_initial(randcov_initial, dedupe_names = TRUE)
  }

  # store partition matrix list
  # this partition_matrix (zeroing covariance across partition_factor groups)
  # is distinct from the "local" big-data partitioning of obdata_list above
  if (!is.null(local$partition_factor)) {
    partition_list <- lapply(obdata_list, function(x) partition_matrix(local$partition_factor, x))
  } else {
    partition_list <- NULL
  }

  # store order
  # records where each original row landed after splitting by local$index, so
  # partition-level results can later be reassembled into original row order
  order <- unlist(split(seq_len(n), local$index), use.names = FALSE)

  # return appropriate list
  list(
    anisotropy = anisotropy, contrasts = dots$contrasts, crs = crs,
    dim_coords = dim_coords, family = family, formula = formula, is_sf = is_sf, local_index = local$index,
    obdata = obdata, obdata_list = obdata_list,
    observed_index = observed_index, offset = offset, ones_list = ones_list, order = order, n = n,
    max_halfdist = max_halfdist, missing_index = missing_index, ncores = local$ncores,
    newdata = newdata, p = p, parallel = local$parallel,
    partition_factor_initial = partition_factor, partition_factor = local$partition_factor,
    partition_list = partition_list, randcov_initial = randcov_initial,
    randcov_list = randcov_list, randcov_names = randcov_names,
    sf_column_name = sf_column_name, size = size, terms = terms_val, var_adjust = local$var_adjust,
    X_list = X_list, xcoord = xcoord, xlevels = xlevels, y_list = y_list, ycoord = ycoord,
    ycoord_orig_name = ycoord_orig_name, ycoord_orig_val = ycoord_orig_val, s2 = s2, diagtol = diagtol,
    range_constrain = range_constrain, range_constrain_value = range_constrain_value
  )
}

#' Build the internal data object used to fit an \code{spgautor()} model
#'
#' @param formula A formula
#' @param family The response family
#' @param data A data frame or \code{sf} object
#' @param spcov_initial A \code{spcov_initial} object
#' @param estmethod The estimation method
#' @param W A neighbor weight matrix (or \code{NULL} to build one internally)
#' @param M A diagonal weighting matrix for the CAR covariance (or \code{NULL})
#' @param random A random effect formula (or \code{NULL})
#' @param randcov_initial A \code{randcov_initial} object (or \code{NULL})
#' @param partition_factor A partition factor formula (or \code{NULL})
#' @param row_st Whether to row-standardize \code{W}
#' @param range_positive Whether the range parameter is constrained positive
#' @param cutoff The neighbor cutoff distance (when \code{W} is built internally)
#' @param ... Additional arguments
#'
#' @return A data object containing the design matrix, response vector,
#'   neighbor structure, and other quantities needed throughout model fitting
#'   and prediction
#'
#' @noRd
get_data_object_spgautor <- function(formula, family, data, spcov_initial,
                                     estmethod, W, M, random, randcov_initial,
                                     partition_factor, row_st, range_positive, cutoff, ...) {
  ## convert sp to sf object
  check_sp_not_supported(data)

  sf_info <- get_areal_sf_info(data)
  is_sf <- sf_info$is_sf
  sf_column_name <- sf_info$sf_column_name
  crs <- sf_info$crs

  # expanding "." in formula
  formula <- expand_formula_dot(formula, data, if (is_sf) sf_column_name else character(0))

  car_neighbor <- build_car_neighbor_structure(data, spcov_initial, W, M, row_st, range_positive, cutoff)
  W <- car_neighbor$W
  W_rowsums <- car_neighbor$W_rowsums
  is_W_connected <- car_neighbor$is_W_connected
  M <- car_neighbor$M
  rho_lb <- car_neighbor$rho_lb
  rho_ub <- car_neighbor$rho_ub

  # subsetting by na and not na values
  ## find response variabale name
  # a cbind(successes, failures) binomial response comes back as a 2-column
  # matrix, so rowSums() collapses it to detect any NA in either column
  response_index <- model.response(model.frame(formula, data, na.action = na.pass))
  if (NCOL(response_index) == 2) {
    na_index <- is.na(rowSums(response_index)) # will be NA if successes or failures NA
  } else {
    na_index <- is.na(response_index)
  }

  # get indices
  observed_index <- which(!na_index)
  missing_index <- which(na_index)
  if (any(na_index)) {
    ## find newdata to be used in prediction later
    newdata <- data[missing_index, , drop = FALSE]
    ## subset original data
    obdata <- data[observed_index, , drop = FALSE]
  } else {
    obdata <- data
    newdata <- NULL
  }


  # finding model frame
  obdata_model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass)
  # finding contrasts as ...
  dots <- list(...)
  if (!"contrasts" %in% names(dots)) {
    dots$contrasts <- NULL
  }
  # model matrix with potential NA
  X <- model.matrix(formula, obdata_model_frame, contrasts = dots$contrasts)
  # finding rows w/out NA
  ob_predictors <- complete.cases(X)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }


  # store X and y
  obdata_model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.omit)
  # store terms
  terms_val <- terms(obdata_model_frame)
  dots <- list(...)
  if (!"contrasts" %in% names(dots)) {
    dots$contrasts <- NULL
  }
  X <- model.matrix(formula, obdata_model_frame, contrasts = dots$contrasts)
  # find induced contrasts and xlevels
  dots$contrasts <- attr(X, "contrasts")
  xlevels <- .getXlevels(terms_val, obdata_model_frame)
  # size is the number of binomial trials per observation; see
  # get_data_object_spglm() above for the full explanation of both response
  # encodings handled here
  y_modr <- model.response(obdata_model_frame)
  if (NCOL(y_modr) == 2) {
    y <- y_modr[, 1, drop = FALSE]
    size <- rowSums(y_modr)
  } else {
    if (family == "binomial") {
      if (is.factor(y_modr)) {
        if (length(levels(y_modr)) != 2) {
          stop("When family is binomial, a factor response must have exactly two levels.", call. = FALSE)
        }
        y_modr <- ifelse(y_modr == levels(y_modr)[1], 0, 1)
      }
      if (is.logical(y_modr)) {
        y_modr <- ifelse(y_modr, 1, 0) # or as.numeric()
      }
      size <- rep(1, NROW(obdata))
    } else {
      size <- NULL
    }
    y <- as.matrix(y_modr, ncol = 1)
  }

  # handle offset
  offset <- model.offset(obdata_model_frame)
  if (!is.null(offset)) {
    offset <- as.matrix(offset, ncol = 1)
  }

  check_response_numeric_and_variable(y)

  # checks on y
  response_checks_glm(family, y, size)

  # store n, p, and ones
  n <- NROW(obdata)
  p <- as.numeric(Matrix::rankMatrix(X, method = "qr"))
  check_rank_collinearity(p, X)
  ones <- matrix(1, nrow = n, ncol = 1)

  check_p_n(p, n, "spautor")

  # find s2 for initial values
  # log(y + 1) approximates the link-scale response; OLS residual variance
  # here is only a starting value for the covariance optimizer, not the fit
  y_trans <- log(y + 1)
  qr_val <- qr(X)
  R_val <- qr.R(qr_val)
  betahat <- backsolve(R_val, qr.qty(qr_val, y_trans))
  resid <- y_trans - X %*% betahat
  s2 <- sum(resid^2) / (n - p)
  diagtol <- 0

  # store random effects list
  if (is.null(random)) {
    randcov_initial <- NULL
    randcov_names <- NULL
    randcov_Zs <- NULL
  } else {
    randcov_names <- get_randcov_names(random)
    randcov_Zs <- get_randcov_Zs(data, randcov_names)
    randcov_initial <- validate_randcov_initial(randcov_initial, dedupe_names = FALSE)
  }

  # coerce to factor
  partition_factor <- coerce_partition_factor(partition_factor, obdata)

  # store partition matrix list
  # zeroes out covariance between observations in different partition_factor
  # groups; layered on top of (not a replacement for) the W neighbor structure
  if (!is.null(partition_factor)) {
    partition_matrix <- partition_matrix(partition_factor, data)
  } else {
    partition_matrix <- NULL
  }

  list(
    anisotropy = FALSE, contrasts = dots$contrasts, crs = crs, family = family,
    formula = formula, data = data, is_sf = is_sf, is_W_connected = is_W_connected,
    missing_index = missing_index, n = n,
    obdata = obdata, observed_index = observed_index, offset = offset, ones = ones, newdata = newdata, p = p,
    partition_factor = partition_factor, partition_matrix = partition_matrix,
    randcov_initial = randcov_initial, randcov_names = randcov_names, randcov_Zs = randcov_Zs,
    sf_column_name = sf_column_name, size = size, terms = terms_val, W = W, W_rowsums = W_rowsums, M = M,
    rho_lb = rho_lb, rho_ub = rho_ub,
    X = X, y = y, xlevels = xlevels, s2 = s2, diagtol = diagtol
  )
}
