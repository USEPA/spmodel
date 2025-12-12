#' Apply the Spatial Decorrelation Transformation to a Newdata Object for Prediction
#'
#' @description Apply the spatial decorrelation transformation to a newdata object.
#'   This object contains explanatory variables that are transformed for prediction
#'   accoring to some spatial decorrelation transformation.
#'
#' @param object A [decorrelate_data()] object.
#' @param newdata A data frame or \code{sf} object in which to
#'   look for variables with which to predict. If a data frame, \code{newdata}
#'   must contain all variables used by \code{formula(object)} and all variables
#'   representing coordinates. If an \code{sf} object, \code{newdata} must contain
#'   all variables used by \code{formula(object)} and coordinates are obtained
#'   from the geometry of \code{newdata}. If omitted, missing data from the
#'   fitted model object are used.
#' @param local A optional logical or list controlling the big data approximation.
#'   If omitted, \code{local} is set
#'   to \code{TRUE} or \code{FALSE} based on the sample size (the number of
#'   non-missing observations in \code{data}) -- if the sample size exceeds 5,000,
#'   \code{local} is set to \code{TRUE}. Otherwise it is set to \code{FALSE}.
#'   If \code{local} is \code{FALSE}, no big data approximation
#'   is implemented. If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{method}: The big data approximation method. If \code{method = "all"},
#'       all observations are used and \code{size} is ignored. If \code{method = "distance"},
#'       the \code{size} data observations closest (in terms of Euclidean distance)
#'       to the observation requiring prediction are used.
#'       If \code{method = "covariance"}, the \code{size} data observations
#'       with the highest covariance with the observation requiring prediction are used.
#'       If random effects and partition factors are not used in estimation and
#'       the spatial covariance function is monotone decreasing,
#'       \code{"distance"} and \code{"covariance"} are equivalent. The default
#'       is \code{"covariance"}.
#'     \item \code{size}: The number of data observations to use when \code{method}
#'       is \code{"distance"} or \code{"covariance"}. The default is 30.
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. This can significantly speed
#'       up computations even when \code{method = "all"} (i.e., no big data
#'       approximation is used), as predictions
#'       are spread out over multiple cores. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 30, method = "covariance", parallel = FALSE)}.
#' @param ... Other arguments.
#'
#' @return A list with many elements that store information about
#'   the fitted model object. Importantly, the list contains the following element:
#'   \itemize{
#'     \item \code{X_newdata}: The original fixed effects design matrix (of explanatory variables) for the prediction data.
#'     \item \code{tX_newdata}: The spatially decorrelated fixed effects design matrix for the prediction data.
#'   }
#'
#' @export
#'
#' @examples
#' params <- spcov_params("exponential", de = 1, ie = 0.2, range = 1e5)
#' decorr <- decorrelate_data(log_cond ~ temp, data = lake, spcov_params = params)
#' decorr_newdata <- decorrelate_newdata(decorr, newdata = lake_preds)
#' head(decorr_newdata$tX_newdata)
decorrelate_newdata <- function(object, newdata, local, ...) {


  if (!inherits(object, "decorrelate_data")) {
    stop("object must have class \"decorrelate_data\".", call. = FALSE)
  }

  if (!missing(local)) {
    object$local <- get_local_list_decorrelate(local)
  }

  # rename relevant quantities
  obdata <- object$obdata
  xcoord <- object$xcoord
  ycoord <- object$ycoord

  if (missing(newdata)) {
    newdata <- object$newdata
  }

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }

  # save spcov param vector
  spcov_params_val <- object$coefficients$spcov

  # save randcov param vector
  randcov_params_val <- object$coefficients$randcov

  # partition factor
  partition_factor_val <- object$partition_factor

  attr_sp <- attr(class(newdata), "package")
  if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
    stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
  }

  if (inherits(newdata, "sf")) {
    newdata <- suppressWarnings(sf::st_centroid(newdata))

    newdata <- sf_to_df(newdata)
    names(newdata)[[which(names(newdata) == ".xcoord")]] <- as.character(xcoord) # only relevant if newdata is sf data is not
    names(newdata)[[which(names(newdata) == ".ycoord")]] <- as.character(ycoord) # only relevant if newdata is sf data is not
  }

  # add back in zero column to cover anisotropy (should make anisotropy only available 1-d)
  if (object$dim_coords == 1) {
    obdata[[ycoord]] <- 0
    newdata[[ycoord]] <- 0
  }

  if (object$anisotropy) { # could just do rotate != 0 || scale != 1
    obdata_aniscoords <- transform_anis(obdata, xcoord, ycoord,
                                        rotate = spcov_params_val[["rotate"]],
                                        scale = spcov_params_val[["scale"]]
    )
    obdata[[xcoord]] <- obdata_aniscoords$xcoord_val
    obdata[[ycoord]] <- obdata_aniscoords$ycoord_val
    object$obdata <- obdata
    newdata_aniscoords <- transform_anis(newdata, xcoord, ycoord,
                                         rotate = spcov_params_val[["rotate"]],
                                         scale = spcov_params_val[["scale"]]
    )
    newdata[[xcoord]] <- newdata_aniscoords$xcoord_val
    newdata[[ycoord]] <- newdata_aniscoords$ycoord_val
  }

  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  if (any(grepl("nmatrix.", attributes(formula_newdata)$dataClasses, fixed = TRUE)) && NROW(newdata) == 1) {
    newdata <- newdata[c(1, 1), , drop = FALSE]
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
    newdata_model <- newdata_model[1, , drop = FALSE]
    # find offset
    offset <- model.offset(newdata_model_frame)
    if (!is.null(offset)) {
      offset <- offset[1]
    }
    newdata <- newdata[1, , drop = FALSE]
  } else {
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    # assumes that predicted observations are not outside the factor levels
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
  }
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(object$X)) # colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # finding rows w/out NA
  ob_predictors <- complete.cases(newdata_model)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }

  # storing newdata as a list
  newdata_rows_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_list <- mapply(x = newdata_rows_list, y = newdata_model_list, FUN = function(x, y) list(row = x, x0 = y), SIMPLIFY = FALSE)

  # randcov stuff
  extra_randcov_list <- get_extra_randcov_list(object, obdata, newdata)
  # reform_bar2_list <- extra_randcov_list$reform_bar2_list
  # Z_index_obdata_list <- extra_randcov_list$Z_index_obdata_list
  # reform_bar1_list <- extra_randcov_list$reform_bar1_list
  # Z_val_obdata_list <- extra_randcov_list$Z_val_obdata_list

  # partition stuff
  extra_partition_list <- get_extra_partition_list(object, obdata, newdata)
  # reform_bar2 <- extra_partition_list$reform_bar2
  # partition_index_obdata <- extra_partition_list$partition_index_obdata

  if (object$local$method == "all") {
    if (object$anisotropy) object$anisotropy <- FALSE # reset anisotropy to
    # FALSE because coordinates already transformed and covmatrix() will rotate/scale them again unnecessarily
    cov_mat <- covmatrix.splm(object)
    cor_mat <- cov_mat / object$total_var
    cor_lowchol <- t(chol(cor_mat))
    rSqrtSigInv_X <- forwardsolve(cor_lowchol, object$X)
    rSqrtSigInv_y <- forwardsolve(cor_lowchol, object$y)
    cor_lowchol_list <- list(
      cor_lowchol = cor_lowchol,
      rSqrtSigInv_X = rSqrtSigInv_X,
      rSqrtSigInv_y = rSqrtSigInv_y
    )
  } else {
    cor_lowchol_list <- NULL
  }

  if (object$local$parallel) {
    cl <- parallel::makeCluster(object$local$ncores)
    output <- parallel::parLapply(cl, newdata_list, get_decorrelate_newdata,
                                  object, cor_lowchol_list, extra_randcov_list, extra_partition_list)
    cl <- parallel::stopCluster(cl)
  } else {
    output <- lapply(newdata_list, get_decorrelate_newdata,
                                  object, cor_lowchol_list, extra_randcov_list, extra_partition_list)
  }

  tX_newdata <- do.call("rbind", lapply(output, function(x) x$tX_newdata))
  rownames(tX_newdata) <- rownames(newdata)
  colnames(tX_newdata) <- colnames(object$X)
  yscale <- do.call("c", lapply(output, function(x) x$yscale))
  names(yscale) <- rownames(newdata)
  yoffset <- do.call("c", lapply(output, function(x) x$yoffset))
  names(yoffset) <- rownames(newdata)
  # remove model matrix structure
  X_newdata <- rbind(newdata_model)
  rownames(X_newdata) <- rownames(newdata)
  colnames(X_newdata) <- colnames(object$X)

  output <- list(
    X_newdata = X_newdata,
    tX_newdata = tX_newdata,
    yscale = yscale,
    yoffset = yoffset
  )
  new_output <- structure(output, class = "decorrelate_newdata")
  new_output

}

get_decorrelate_newdata <- function(newdata_list, object, cor_lowchol_list, extra_randcov_list, extra_partition_list) {

  obdata <- object$obdata
  xcoord <- object$xcoord
  ycoord <- object$ycoord
  X <- object$X
  y <- object$y

  # storing partition vector
  if (!is.null(object$partition_factor)) {
    partition_vector <- partition_vector(object$partition_factor,
                                         data = object$obdata,
                                         newdata = newdata_list$row, reform_bar2 = extra_partition_list$reform_bar2,
                                         partition_index_data = extra_partition_list$partition_index_obdata
    )
  } else {
    partition_vector <- NULL
  }

  dist_vector <- spdist_vectors(newdata_list$row, obdata, xcoord, ycoord, object$dim_coords)

  # making random vector if necessary
  if (!is.null(object$random)) {
    randcov_vector_val <- randcov_vector(object$coefficients$randcov, object$obdata, newdata_list$row,
                                         extra_randcov_list$reform_bar2_list, extra_randcov_list$Z_index_obdata_list)
  } else {
    randcov_vector_val <- NULL
  }

  # making the covariance vector
  cov_vector_val <- cov_vector(object$coefficients$spcov, dist_vector, randcov_vector_val, partition_vector)
  cov_vector_val <- as.numeric(cov_vector_val)

  # subsetting data if method distance
  if (object$local$method == "distance") {
    n <- length(cov_vector_val)
    # want the smallest distance here and order goes from smallest first to largest last (keep last values with are smallest distance)
    nn_index <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, object$local$size))]
    obdata <- obdata[nn_index, , drop = FALSE]
    X <- X[nn_index, , drop = FALSE]
    y <- y[nn_index]
    cov_vector_val <- cov_vector_val[nn_index]
  }

  if (object$local$method == "covariance") {
    n <- length(cov_vector_val)
    # want the largest covariance here and order goes from smallest first to largest last (keep last values which are largest covariance)
    cov_index <- order(as.numeric(cov_vector_val))[seq(from = n, to = max(1, n - object$local$size + 1))] # use abs() here?
    obdata <- obdata[cov_index, , drop = FALSE]
    X <- X[cov_index, , drop = FALSE]
    y <- y[cov_index]
    cov_vector_val <- cov_vector_val[cov_index]
  }

  if (object$local$method %in% c("distance", "covariance")) {
    if (!is.null(object$random)) {
      randcov_names <- get_randcov_names(object$random)
      xlev_list <- lapply(extra_randcov_list$Z_index_obdata_list, function(x) x$reform_bar2_xlev)
      randcov_Zs <- get_randcov_Zs(obdata, randcov_names, xlev_list = xlev_list)
    }
    partition_matrix_val <- partition_matrix(object$partition_factor, obdata)
    cov_matrix_val <- cov_matrix(
      object$coefficients$spcov, spdist(obdata, xcoord, ycoord), object$coefficients$randcov,
      randcov_Zs, partition_matrix_val,
      diagtol = object$diagtol
    )
    cor_matrix_val <- cov_matrix_val / object$total_var
    cor_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(cor_matrix_val)))
  } else {
    cor_lowchol <- cor_lowchol_list$cor_lowchol
  }

  cor_vector_val <- cov_vector_val / object$total_var

  rSqrtSigInv_r0 <- forwardsolve(cor_lowchol, cor_vector_val)
  r0_SigInv_r0 <- crossprod(rSqrtSigInv_r0, rSqrtSigInv_r0)
  w <- as.numeric(1 - r0_SigInv_r0)

  if (object$local$method %in% c("distance", "covariance")) {
    rSqrtSigInv_X <- forwardsolve(cor_lowchol, X)
    rSqrtSigInv_y <- forwardsolve(cor_lowchol, y)
  } else {
    rSqrtSigInv_X <- cor_lowchol_list$rSqrtSigInv_X
    rSqrtSigInv_y <- cor_lowchol_list$rSqrtSigInv_y
  }


  sqrt_w <- sqrt(w)
  tX_newdata <- (newdata_list$x0 - crossprod(rSqrtSigInv_r0, rSqrtSigInv_X)) / sqrt_w
  yoffset <- crossprod(rSqrtSigInv_r0, rSqrtSigInv_y)

  list(
    tX_newdata = tX_newdata,
    yscale = sqrt_w,
    yoffset = as.numeric(yoffset)
  )


}
