#' Create a covariance matrix
#'
#' Create a covariance matrix from a fitted model object.
#'
#' @param object A fitted model object (e.g., [splm()], [spautor()], [spglm()], or [spgautor()]).
#' @param newdata If omitted, the covariance matrix of
#'   the observed data is returned. If provided, \code{newdata} is
#'   a data frame or \code{sf} object that contains coordinate information
#'   required to construct the covariance between \code{newdata} and
#'   the observed data. If a data frame, \code{newdata}
#'   must contain variables that represent coordinates having the same name as
#'   the coordinates from the observed data used to fit \code{object}. If an
#'   \code{sf} object, coordinates are obtained from the geometry of \code{newdata}.
#' @param cov_type The type of covariance matrix returned. If \code{newdata}
#'   is omitted or \code{cov_type} is \code{"obs.obs"},
#'   the \eqn{n \times n} covariance matrix of the observed
#'   data is returned, where \eqn{n} is the sample size used to fit \code{object}.
#'   If \code{newdata} is provided and \code{cov_type} is \code{"pred.obs"}
#'   (the default when \code{newdata} is provided),
#'   the \eqn{m \times n} covariance matrix of the prediction and observed data is returned,
#'   where \eqn{m} is the number of elements in the prediction data.
#'   If \code{newdata} is provided and \code{cov_type} is \code{"obs.pred"},
#'   the \eqn{n \times m} covariance matrix of the observed and prediction data is returned.
#'   If \code{newdata} is provided and \code{cov_type} is \code{"pred.pred"},
#'   the \eqn{m \times m} covariance matrix of the prediction data is returned.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return If \code{newdata} is omitted, the covariance matrix of the observed
#'   data, which has dimension n x n, where n is the sample size used to fit \code{object}.
#'   If \code{newdata} is provided, the covariance matrix between the unobserved (new)
#'   data and the observed data, which has dimension m x n, where m is the number of
#'   new observations and n is the sample size used to fit \code{object}.
#'
#' @order 1
#' @export
#'
#' @examples
#' spmod <- splm(z ~ water + tarp,
#'   data = caribou,
#'   spcov_type = "exponential", xcoord = x, ycoord = y
#' )
#' covmatrix(spmod)
covmatrix <- function(object, newdata, ...) {
  UseMethod("covmatrix", object)
}
#' @rdname covmatrix
#' @method covmatrix splm
#' @order 2
#' @export
covmatrix.splm <- function(object, newdata, cov_type, ...) {

  if (missing(newdata)) {
    cov_type <- "obs.obs"
  } else {
    if (missing(cov_type)) {
      cov_type <- "pred.obs"
    }
  }

  if (cov_type != "obs.obs") {
    if (is.null(newdata)) {
      stop("No prediction data for which to create a covariance matrix.", call. = FALSE)
    }
  }

  if (cov_type == "pred.pred") {
    if (inherits(object$newdata, "sf")) {
      object$obdata <- sf_to_df(object$newdata)
    } else {
      object$obdata <- newdata
    }
    return(covmatrix(object))
  }

  # spcov
  spcov_params_val <- coef(object, type = "spcov")

  # randcov
  randcov_params_val <- coef(object, type = "randcov")

  # if (missing(newdata)) {
  if (cov_type == "obs.obs") {
    # coordinate stuff
    if (object$anisotropy) {
      new_coords <- transform_anis(
        object$obdata, object$xcoord, object$ycoord,
        spcov_params_val[["rotate"]], spcov_params_val[["scale"]]
      )
      dist_matrix <- spdist(xcoord_val = new_coords$xcoord_val, ycoord_val = new_coords$ycoord_val)
    } else {
      dist_matrix <- spdist(object$obdata, object$xcoord, object$ycoord)
    }

    # random effects
    randcov_names <- get_randcov_names(object$random)
    randcov_Zs_val <- get_randcov_Zs(object$obdata, randcov_names)
    randcov_matrix_val <- randcov_matrix(randcov_params_val, randcov_Zs = randcov_Zs_val)

    # partition factor
    partition_matrix_val <- partition_matrix(object$partition_factor, object$obdata)

    # cov matrix
    cov_val <- cov_matrix(spcov_params_val, dist_matrix, randcov_params_val,
      randcov_Zs_val,
      partition_matrix = partition_matrix_val,
      diagtol = object$diagtol
    )
  } else if (cov_type %in% c("pred.obs", "obs.pred")) {

    # rename relevant quantities
    obdata <- object$obdata
    xcoord <- object$xcoord
    ycoord <- object$ycoord

    # transform newdata if required
    if (inherits(newdata, "sf")) {
      newdata <- suppressWarnings(sf::st_centroid(newdata))

      newdata <- sf_to_df(newdata)
      names(newdata)[[which(names(newdata) == ".xcoord")]] <- as.character(xcoord) # only relevant if newdata is sf data is not
      names(newdata)[[which(names(newdata) == ".ycoord")]] <- as.character(ycoord) # only relevant if newdata is sf data is not
    }

    # transform aniosotropy coordinates
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
      newdata_aniscoords <- transform_anis(newdata, xcoord, ycoord,
        rotate = spcov_params_val[["rotate"]],
        scale = spcov_params_val[["scale"]]
      )
      newdata[[xcoord]] <- newdata_aniscoords$xcoord_val
      newdata[[ycoord]] <- newdata_aniscoords$ycoord_val
    }

    # coordinate stuff
    dist_vector <- spdist_vectors(newdata, obdata, xcoord, ycoord, object$dim_coords)

    # random effects
    randcov_vector_val <- randcov_vector(randcov_params_val, obdata, newdata)

    # partition factor
    partition_vector_val <- partition_vector(object$partition_factor, obdata, newdata)

    # cov vector
    cov_val <- cov_vector(spcov_params_val, dist_vector, randcov_vector_val, partition_vector_val)

    if (cov_type == "obs.pred") {
      cov_val <- t(cov_val)
    }

  } else {
    stop('cov_type must be "obs.obs", "obs.pred", "pred.obs", "pred.pred"', call. = FALSE)
  }

  # return covariance value as a base R matrix (not a Matrix matrix)
  as.matrix(cov_val)
}

#' @rdname covmatrix
#' @method covmatrix spautor
#' @order 3
#' @export
covmatrix.spautor <- function(object, newdata, cov_type, ...) {

  if (missing(newdata)) {
    cov_type <- "obs.obs"
  } else {
    if (missing(cov_type)) {
      cov_type <- "pred.obs"
    }
  }

  if (cov_type != "obs.obs") {
    if (is.null(newdata)) {
      stop("No prediction data for which to create a covariance matrix.", call. = FALSE)
    }
  }

  # spcov
  spcov_params_val <- coef(object, type = "spcov")

  # randcov
  randcov_params_val <- coef(object, type = "randcov")

  # weights matrix
  dist_matrix <- object$W

  # random effect
  randcov_names <- get_randcov_names(object$random)
  randcov_Zs_val <- get_randcov_Zs(object$data, randcov_names)
  randcov_matrix_val <- randcov_matrix(randcov_params_val, randcov_Zs = randcov_Zs_val)

  # partition factor
  partition_matrix_val <- partition_matrix(object$partition_factor, object$data)

  # full cov matrix
  cov_matrix_val <- cov_matrix(spcov_params_val, dist_matrix, randcov_params_val,
    randcov_Zs_val,
    partition_matrix = partition_matrix_val, object$M
  )

  if (cov_type == "obs.obs") {
    # observed covariance matrix
    cov_val <- cov_matrix_val[object$observed_index, object$observed_index, drop = FALSE]
  } else if (cov_type == "pred.obs") {
    # new covariance vector
    cov_val <- cov_matrix_val[object$missing_index, object$observed_index, drop = FALSE]
  } else if (cov_type == "obs.pred") {
    # new covariance vector
    cov_val <- cov_matrix_val[object$observed_index, object$missing_index, drop = FALSE]
  } else if (cov_type == "pred.pred") {
    cov_val <- cov_matrix_val[object$missing_index, object$missing_index, drop = FALSE]
  } else {
    stop('cov_type must be "obs.obs", "obs.pred", "pred.obs", "pred.pred"', call. = FALSE)
  }

  # return covariance value as a base R matrix (not a Matrix matrix)
  as.matrix(cov_val)
}
