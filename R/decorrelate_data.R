#' Apply the Spatial Decorrelation Transformation to a Data Object
#'
#' @description Apply the spatial decorrelation transformation to a data object.
#'   This object contains the transformed explanatory and response variables
#'   which can be used to fit a machine learning model. This object also contains
#'   information needed to decorrelate prediction data.
#'
#' @inheritParams decorrelate
#'
#' @details The spatial decorrelation transformation is a preprocessing transformation
#'   that reduces the impacts of spatial dependence (i.e., covariance, correlation)
#'   on machine learning models. See [decorrelate()] and Heaton et al., 2025 for more details.
#'
#'
#' @return A list with many elements that store information about
#'   the fitted model object. Importantly, the list contains the following elements:
#'   \itemize{
#'     \item \code{X}: The original fixed effects design matrix (of explanatory variables)
#'     \item \code{y}: The original response variable
#'     \item \code{tX}: The spatially decorrelated transformed fixed effects design matrix
#'     \item \code{ty}: The spatially decorrelated transformed response variable
#'   }
#'
#' @export
#'
#' @seealso [decorrelate()] [spcov_params()] [randcov_params()]
#'
#'
#' @references Matthew J. Heaton, Andrew Millane, and Jake S. Rhodes. 2025. A Scalable
#'   Spatial Decorrelation Preprocessing Approach for Machine and Deep Learning.
#'   \emph{Journal of Data Science}. 1-15, DOI 10.6339/25-JDS1210
#'
#' @examples
#' params <- spcov_params("exponential", de = 1, ie = 0.2, range = 1e5)
#' decorr <- decorrelate_data(log_cond ~ temp, data = lake, spcov_params = params)
#' head(cbind(decorr$X, decorr$tX))
#' head(cbind(decorr$y, decorr$ty))
decorrelate_data <- function(formula, data, spcov_params, xcoord, ycoord, randcov_params, partition_factor, ordering, local, ...) {

  # non standard evaluation for x and y coordinates (only meaningful at this,
  # the direct calling frame -- decorrelate_data_internal() receives the
  # already-substituted value and must not re-substitute)
  xcoord <- substitute(xcoord)
  ycoord <- substitute(ycoord)

  # set randcov_initial NULL if necessary. Unlike decorrelate_data_internal()
  # (which always receives an already-resolved `random` from its caller),
  # decorrelate_data() has no `random` argument of its own and must derive it
  # from randcov_params.
  if (missing(randcov_params) || is.null(randcov_params)) {
    random <- NULL
    randcov_params <- NULL
  } else {
    random <- reformulate(names(randcov_params))
  }

  # decorrelate_data_internal() uses is.null(ordering)/is.null(local), which
  # forces evaluation, so a truly-missing argument must be resolved to an
  # explicit NULL here rather than passed through by bare symbol.
  if (missing(partition_factor)) partition_factor <- NULL
  if (missing(ordering)) ordering <- NULL
  if (missing(local)) local <- NULL

  decorrelate_data_internal(
    formula = formula,
    data = data,
    spcov_params = spcov_params,
    xcoord = xcoord,
    ycoord = ycoord,
    random = random,
    randcov_params = randcov_params,
    partition_factor = partition_factor,
    ordering = ordering,
    local = local,
    ...
  )
}

decorrelate_data_internal_part1 <- function(formula, data, spcov_type, xcoord, ycoord, random, partition_factor, ordering, local, ...) {

  if (missing(random)) {
    random <- NULL
  }

  # set partition factor if necessary
  if (missing(partition_factor) || is.null(partition_factor)) {
    partition_factor <- NULL
  }

  # non standard evaluation for x and y coordinates
  # xcoord <- substitute(xcoord)
  # ycoord <- substitute(ycoord)

  # get data object
  data_object <- get_data_object_splm(
    formula = formula,
    data = data,
    spcov_initial = spcov_initial(spcov_type), # default placeholder
    xcoord = xcoord,
    ycoord = ycoord,
    estmethod = "reml",  # default placeholder
    anisotropy = FALSE, # default placeholder
    random = random,
    randcov_initial = NULL, # default placeholder
    partition_factor = NULL, # default placeholder
    local = FALSE, # default placeholder
    range_constrain = FALSE, # default placeholder
    ...
  )

  if (missing(local)) {
    local <- NULL
  }
  if (is.null(local)) {
    if (data_object$n > 5000) {
      local <- TRUE
      message("Because the sample size exceeds 5,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }
  local <- get_local_list_decorrelate(local)

  if (is.null(partition_factor)) {
    partition_matrix_val <- NULL
  } else {
    partition_matrix_val <- partition_matrix(partition_factor, data)
  }

  obdata <- data_object$obdata
  X <- data_object$X_list[[1]]
  y <- data_object$y_list[[1]]
  xcoord_val <- obdata[[data_object$xcoord]]
  ycoord_val <- obdata[[data_object$ycoord]]

  index <- seq(1, data_object$n)

  # do ordering here
  if (is.null(ordering)) {
    any_dup <- any(duplicated(cbind(xcoord_val, ycoord_val)))
    if (any_dup) {
      ordering <- "grts"
    } else {
      ordering <- "maxmin"
    }
  }

  if (!ordering %in% c("middleout", "outsidein", "coordinate", "maxmin", "grts", "random", "none")) {
    stop("Invalid ordering argument. Argument must be \"middleout\", \"outsidein\", \"coordinate\", \"maxmin\", \"grts\", \"random\", or \"none\".", call. = FALSE)
  }

  # ordering done with separate random elements (grts, random) for each grid item (fix)
  ord <- get_decorrelate_order(ordering, xcoord_val, ycoord_val)

  # order all values
  X <- X[ord$order, , drop = FALSE]
  y <- y[ord$order, , drop = FALSE]
  xcoord_val <- xcoord_val[ord$order]
  ycoord_val <- ycoord_val[ord$order]
  if (!is.null(partition_matrix_val)) {
    partition_matrix_val <- partition_matrix_val[ord$order, ord$order, drop = FALSE]
  }

  list(
    index = index,
    X = X,
    y = y,
    xcoord_val = xcoord_val,
    ycoord_val = ycoord_val,
    local = local,
    random = random,
    partition_factor = partition_factor,
    partition_matrix_val = partition_matrix_val,
    ord = ord,
    ordering = ordering,
    data_object = data_object
  )
}

decorrelate_data_internal_part2 <- function(spcov_params, randcov_params, decorrelate_part1_object, ...) {

  index <- decorrelate_part1_object$index
  total_var <- decorrelate_part1_object$total_var
  X <- decorrelate_part1_object$X
  y <- decorrelate_part1_object$y
  xcoord_val <- decorrelate_part1_object$xcoord_val
  ycoord_val <- decorrelate_part1_object$ycoord_val
  local <- decorrelate_part1_object$local
  random <- decorrelate_part1_object$random
  partition_factor <- decorrelate_part1_object$partition_factor
  partition_matrix_val <- decorrelate_part1_object$partition_matrix_val
  ord <- decorrelate_part1_object$ord
  ordering <- decorrelate_part1_object$ordering
  data_object <- decorrelate_part1_object$data_object


    if (spcov_params[["rotate"]] != 0 || spcov_params[["scale"]] != 1) {
      anisotropy <- TRUE
    } else {
      anisotropy <- FALSE
    }

    if (anisotropy) {
      obdata_aniscoords <- transform_anis2(xcoord_val, ycoord_val, spcov_params[["rotate"]], spcov_params[["scale"]])
      xcoord_val <- obdata_aniscoords$xcoord
      ycoord_val <- obdata_aniscoords$ycoord
    }

  if (!is.null(randcov_params)) {
    randcov_matrix_val <- randcov_matrix(randcov_params, data_object$randcov_list[[1]])
    randcov_matrix_val <- randcov_matrix_val[ord$order, ord$order, drop = FALSE]
  } else {
    randcov_matrix_val <- NULL
  }


  total_var <- sum(spcov_params[["de"]], spcov_params[["ie"]], randcov_params)

  if (local$parallel) {
    cl <- parallel::makeCluster(local$ncores)
    vals <- parallel::parLapply(cl, index, get_decorrelated_value, spcov_params, total_var, X, y, xcoord_val, ycoord_val, local, randcov_matrix_val, partition_matrix_val)
    cl <- parallel::stopCluster(cl)
  } else {
    vals <- lapply(index, get_decorrelated_value, spcov_params, total_var, X, y, xcoord_val, ycoord_val, local, randcov_matrix_val, partition_matrix_val)
  }
  X <- do.call("rbind", lapply(vals, function(x) x$X))
  y <- do.call("rbind", lapply(vals, function(x) x$y))
  tX <- do.call("rbind", lapply(vals, function(x) x$tX))
  ty <- do.call("rbind", lapply(vals, function(x) x$ty))

  # undo vecchia ordering here
  X <- X[ord$inv_order, , drop = FALSE]
  y <- y[ord$inv_order, , drop = FALSE]
  tX <- tX[ord$inv_order, , drop = FALSE]
  ty <- ty[ord$inv_order, , drop = FALSE]

  coefs <- list(spcov = spcov_params, randcov = randcov_params)
  output <- list(
    obdata = data_object$obdata,
    coefficients = coefs,
    X = X,
    y = as.vector(y),
    tX = tX,
    ty = as.vector(ty),
    xcoord = data_object$xcoord,
    ycoord = data_object$ycoord,
    random = random,
    partition_factor = partition_factor,
    dim_coords = data_object$dim_coords,
    terms = data_object$terms,
    xlevels = data_object$xlevels,
    contrasts = data_object$contrasts,
    local = local,
    newdata = data_object$newdata,
    anisotropy = anisotropy,
    diagtol = data_object$diagtol,
    total_var = total_var,
    ordering = ordering
  )
  new_output <- structure(output, class = "decorrelate_data")
  new_output
}

# non standard evaluation for x and y coordinates
# substitute only works when the function is the parent function, so commenting
# out for use with decorrelate()
decorrelate_data_internal <- function(formula, data, spcov_params, xcoord, ycoord, random, randcov_params, partition_factor, ordering, local, ...) {

  if (spcov_params[["rotate"]] != 0 || spcov_params[["scale"]] != 1) {
    anisotropy <- TRUE
  } else {
    anisotropy <- FALSE
  }

  # set randcov_initial NULL if necessary
  if (missing(randcov_params) || is.null(randcov_params)) {
    random <- NULL
    randcov_params <- NULL
  }

  # set partition factor if necessary
  if (missing(partition_factor) || is.null(partition_factor)) {
    partition_factor <- NULL
  }

  # non standard evaluation for x and y coordinates
  # xcoord <- substitute(xcoord)
  # ycoord <- substitute(ycoord)

  # get data object
  data_object <- get_data_object_splm(
    formula = formula,
    data = data,
    spcov_initial = spcov_initial(class(spcov_params)), # default placeholder
    xcoord = xcoord,
    ycoord = ycoord,
    estmethod = "reml",  # default placeholder
    anisotropy = anisotropy,
    random = random,
    randcov_initial = NULL, # default placeholder
    partition_factor = NULL, # default placeholder
    local = FALSE, # default placeholder
    range_constrain = FALSE, # default placeholder
    ...
  )

  if (missing(local)) {
    local <- NULL
  }
  if (is.null(local)) {
    if (data_object$n > 5000) {
      local <- TRUE
      message("Because the sample size exceeds 5,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }
  local <- get_local_list_decorrelate(local)

  # random effects
  if (is.null(random)) {
    randcov_matrix_val <- NULL
  } else {
    names(randcov_params) <- get_randcov_names(random) # fixes names
    randcov_matrix_val <- randcov_matrix(randcov_params, data_object$randcov_list[[1]])
  }

  if (is.null(partition_factor)) {
    partition_matrix_val <- NULL
  } else {
    partition_matrix_val <- partition_matrix(partition_factor, data)
  }


  # return original information here if spcov_type "none"


  obdata <- data_object$obdata
  X <- data_object$X_list[[1]]
  y <- data_object$y_list[[1]]
  xcoord_val <- obdata[[data_object$xcoord]]
  ycoord_val <- obdata[[data_object$ycoord]]
  if (anisotropy) {
    obdata_aniscoords <- transform_anis2(xcoord_val, ycoord_val, spcov_params[["rotate"]], spcov_params[["scale"]])
    xcoord_val <- obdata_aniscoords$xcoord
    ycoord_val <- obdata_aniscoords$ycoord
  }


  index <- seq(1, data_object$n)
  total_var <- sum(spcov_params[["de"]], spcov_params[["ie"]], randcov_params)

  # do ordering here
  if (is.null(ordering)) {
    any_dup <- any(duplicated(cbind(xcoord_val, ycoord_val)))
    if (any_dup) {
      ordering <- "grts"
    } else {
      ordering <- "maxmin"
    }
  }

  if (!ordering %in% c("middleout", "outsidein", "coordinate", "maxmin", "grts", "random", "none")) {
    stop("Invalid ordering argument. Argument must be \"middleout\", \"outsidein\", \"coordinate\", \"maxmin\", \"grts\", \"random\", or \"none\".", call. = FALSE)
  }

  # ordering done with separate random elements (grts, random) for each grid item (fix)
  ord <- get_decorrelate_order(ordering, xcoord_val, ycoord_val) # remove this eventually replace with ordering_list

  # order all values
  X <- X[ord$order, , drop = FALSE]
  y <- y[ord$order, , drop = FALSE]
  xcoord_val <- xcoord_val[ord$order]
  ycoord_val <- ycoord_val[ord$order]
  if (!is.null(randcov_matrix_val)) {
    randcov_matrix_val <- randcov_matrix_val[ord$order, ord$order, drop = FALSE]
  }
  if (!is.null(partition_matrix_val)) {
    partition_matrix_val <- partition_matrix_val[ord$order, ord$order, drop = FALSE]
  }


  if (local$parallel) {
    cl <- parallel::makeCluster(local$ncores)
    vals <- parallel::parLapply(cl, index, get_decorrelated_value, spcov_params, total_var, X, y, xcoord_val, ycoord_val, local, randcov_matrix_val, partition_matrix_val)
    cl <- parallel::stopCluster(cl)
  } else {
    vals <- lapply(index, get_decorrelated_value, spcov_params, total_var, X, y, xcoord_val, ycoord_val, local, randcov_matrix_val, partition_matrix_val)
  }
  X <- do.call("rbind", lapply(vals, function(x) x$X))
  y <- do.call("rbind", lapply(vals, function(x) x$y))
  tX <- do.call("rbind", lapply(vals, function(x) x$tX))
  ty <- do.call("rbind", lapply(vals, function(x) x$ty))

  # undo vecchia ordering here
  X <- X[ord$inv_order, , drop = FALSE]
  y <- y[ord$inv_order, , drop = FALSE]
  tX <- tX[ord$inv_order, , drop = FALSE]
  ty <- ty[ord$inv_order, , drop = FALSE]

  coefs <- list(spcov = spcov_params, randcov = randcov_params)
  output <- list(
    obdata = data_object$obdata,
    coefficients = coefs,
    X = X,
    y = as.vector(y),
    tX = tX,
    ty = as.vector(ty),
    xcoord = data_object$xcoord,
    ycoord = data_object$ycoord,
    random = random,
    partition_factor = partition_factor,
    dim_coords = data_object$dim_coords,
    terms = data_object$terms,
    xlevels = data_object$xlevels,
    contrasts = data_object$contrasts,
    local = local,
    newdata = data_object$newdata,
    anisotropy = data_object$anisotropy,
    diagtol = data_object$diagtol,
    total_var = total_var,
    ordering = ordering
  )
  new_output <- structure(output, class = "decorrelate_data")
  new_output
}


get_decorrelated_value <- function(index, spcov_params, total_var, X, y, xcoord_val, ycoord_val, local, randcov_matrix, partition_matrix) {

  if (index == 1) {
    X <- X[index, , drop = FALSE]
    y <- y[index, , drop = FALSE]
    return(list(X = X, y = y, tX = X, ty = y))
  }

  index_new <- index
  index_old <- seq(1, index - 1)

  X_new <- X[index_new, , drop = FALSE]
  X_old <- X[index_old, , drop = FALSE]
  y_new <- y[index_new, , drop = FALSE]
  y_old <- y[index_old, , drop = FALSE]

  xcoord_val_new <- xcoord_val[index_new]
  ycoord_val_new <- ycoord_val[index_new]
  xcoord_val_old <- xcoord_val[index_old]
  ycoord_val_old <- ycoord_val[index_old]

  dists_new <- spdist_vectors2(xcoord_val_new, ycoord_val_new, xcoord_val_old, ycoord_val_old)
  # random effects if necessary
  if (is.null(randcov_matrix)) {
    randcov_vector <- NULL
  } else {
    randcov_vector <- randcov_matrix[index_new, index_old, drop = FALSE]
  }
  # partition factor if necessary
  if (is.null(partition_matrix)) {
    partition_vector <- NULL
  } else {
    partition_vector <- partition_matrix[index_new, index_old, drop = FALSE]
  }
  cov_vec_new <- cov_vector(spcov_params, dists_new, randcov_vector = randcov_vector, partition_vector = partition_vector)
  cov_vec_new <- as.numeric(cov_vec_new)

  # random effects if necessary
  if (is.null(randcov_matrix)) {
    randcov_matrix <- NULL
  } else {
    randcov_matrix <- randcov_matrix[index_old, index_old, drop = FALSE]
  }
  # partition factor if necessary
  if (is.null(partition_matrix)) {
    partition_matrix <- NULL
  } else {
    partition_matrix <- partition_matrix[index_old, index_old, drop = FALSE]
  }

  if (local$method == "all" ) { # || index <= local$size
    dists_old <- spdist(xcoord_val = xcoord_val_old, ycoord_val = ycoord_val_old)
  }

  # subsetting data if method distance
  if (local$method == "distance") { # && ndex > local$size
    n <- length(cov_vec_new)
    # want the smallest distance here and order goes from smallest first to largest last (keep last values with are smallest distance)
    nn_index <- order(as.numeric(dists_new))[seq(from = 1, to = min(n, local$size))]
    X_old <- X_old[nn_index, , drop = FALSE]
    y_old <- y_old[nn_index, , drop = FALSE]
    xcoord_val_old <- xcoord_val[nn_index]
    ycoord_val_old <- ycoord_val[nn_index]
    dists_old <- spdist(xcoord_val = xcoord_val_old, ycoord_val = ycoord_val_old)
    if (!is.null(randcov_matrix)) {
      randcov_matrix <- randcov_matrix[nn_index, nn_index, drop = FALSE]
    }
    if (!is.null(partition_matrix)) {
      partition_matrix <- partition_matrix[nn_index, nn_index, drop = FALSE]
    }
    cov_vec_new <- cov_vec_new[nn_index]
  }

  if (local$method == "covariance") { # && index > local$size
    n <- length(cov_vec_new)
    # want the largest covariance here and order goes from smallest first to largest last (keep last values which are largest covariance)
    # TODO: ranks by raw covariance, not |covariance|. For monotone decreasing spcov_types
    # (exponential, spherical, etc.) this is equivalent to nearest-neighbor selection, but for
    # spcov_types with negative lobes (wave, cosine, jbessel) it can pass over strongly
    # negatively-correlated neighbors in favor of weakly positive ones. Revisit whether to rank
    # by abs(cov_vec_new) instead.
    cov_index <- order(as.numeric(cov_vec_new))[seq(from = n, to = max(1, n - local$size + 1))]
    X_old <- X_old[cov_index, , drop = FALSE]
    y_old <- y_old[cov_index, , drop = FALSE]
    xcoord_val_old <- xcoord_val[cov_index]
    ycoord_val_old <- ycoord_val[cov_index]
    dists_old <- spdist(xcoord_val = xcoord_val_old, ycoord_val = ycoord_val_old)
    if (!is.null(randcov_matrix)) {
      randcov_matrix <- randcov_matrix[cov_index, cov_index, drop = FALSE]
    }
    if (!is.null(partition_matrix)) {
      partition_matrix <- partition_matrix[cov_index, cov_index, drop = FALSE]
    }
    cov_vec_new <- cov_vec_new[cov_index]
  }

  cov_mat_old <- cov_matrix2(spcov_params, dist_matrix = dists_old, randcov_matrix = randcov_matrix, partition_matrix = partition_matrix)

  cor_vec_new <- cov_vec_new / total_var
  cor_mat_old <- cov_mat_old / total_var

  rSig_upchol <- Matrix::chol(Matrix::forceSymmetric(cor_mat_old))
  rSig_lowchol <- t(rSig_upchol)

  rSqrtSigInv_r0 <- forwardsolve(rSig_lowchol, cor_vec_new)
  r0_SigInv_r0 <- crossprod(rSqrtSigInv_r0, rSqrtSigInv_r0)
  w <- as.numeric(1 - r0_SigInv_r0)

  rSqrtSigInv_Xold <- forwardsolve(rSig_lowchol, X_old)
  rSqrtSigInv_yold <- forwardsolve(rSig_lowchol, y_old)

  tX_new <- (X_new - crossprod(rSqrtSigInv_r0, rSqrtSigInv_Xold)) / sqrt(w)
  ty_new <- (y_new - crossprod(rSqrtSigInv_r0, rSqrtSigInv_yold)) / sqrt(w)

  list(X = X_new, y = y_new, tX = tX_new, ty = ty_new)
}

get_local_list_decorrelate <- function(local) {


  if (is.logical(local)) {
    if (local) {
      local <- list(method = "covariance", size = 30)
    } else {
      local <- list(method = "all")
    }
  }

  names_local <- names(local)

  # errors
  if ("method" %in% names_local) {
    if (!local$method %in% c("all", "covariance", "distance")) {
      stop("Invalid local method. Local method must be \"all\", \"covariance\", or \"distance\".", call. = FALSE)
    }
  }

  if (!"method" %in% names_local) {
    local$method <- "covariance"
  }

  if (local$method %in% c("distance", "covariance") && !"size" %in% names_local) {
    local$size <- 30
  }

  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$parallel) {
    if (!"ncores" %in% names_local) {
      local$ncores <- parallel::detectCores()
    }
  }

  local

}

get_decorrelate_order <- function(ordering, xcoord_val, ycoord_val) {

  n <- length(xcoord_val)

  if (ordering == "none") {
    ord <- seq(1, n)
  }

  if (ordering == "random") {
    ord <- seq(1, n)
    ord <- sample(ord)
  }

  if (ordering %in% "maxmin") {
    if (!requireNamespace("GPvecchia", quietly = TRUE)) {
      stop("Install the GPvecchia package before using \"maxmin\" ordering", call. = FALSE)
    } else {
      ord <- GPvecchia::order_maxmin_exact(cbind(xcoord_val, ycoord_val))
    }
  }

  if (ordering %in% "middleout") {
    if (!requireNamespace("GPvecchia", quietly = TRUE)) {
      stop("Install the GPvecchia package before using \"middleout\" ordering", call. = FALSE)
    } else {
      ord <- GPvecchia::order_middleout(cbind(xcoord_val, ycoord_val))
    }
  }

  if (ordering %in% "outsidein") {
    if (!requireNamespace("GPvecchia", quietly = TRUE)) {
      stop("Install the GPvecchia package before using \"outsidein\" ordering", call. = FALSE)
    } else {
      ord <- GPvecchia::order_outsidein(cbind(xcoord_val, ycoord_val))
    }
  }

  if (ordering %in% "coordinate") {
    if (!requireNamespace("GPvecchia", quietly = TRUE)) {
      stop("Install the GPvecchia package before using \"coordinate\" ordering", call. = FALSE)
    } else {
      ord <- GPvecchia::order_coordinate(cbind(xcoord_val, ycoord_val), coordinate = c(1, 2))
    }
  }

  if (ordering == "grts") {
    if (!requireNamespace("spsurvey", quietly = TRUE)) {
      stop("Install the spsurvey package before using \"grts\" ordering", call. = FALSE)
    } else {
      dat <- data.frame(xcoord_val = xcoord_val, ycoord_val = ycoord_val, ord = seq(1, n))
      sframe <- st_as_sf(dat, coords = c("xcoord_val", "ycoord_val"), crs = NA)
      samp <- spsurvey::grts(sframe, n_base = n, projcrs_check = FALSE)
      ord <- samp$sites_base$ord
    }
  }

  list(order = ord, inv_order = order(ord))
}
