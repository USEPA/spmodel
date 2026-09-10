#' Build the prediction setup object used by \code{predict.splm()}
#'
#' @param object A fitted model object from [splm()]
#' @param newdata Data requiring prediction (possibly missing, to predict \code{object$newdata})
#' @param scale A scale multiplier for prediction standard errors (or \code{NULL})
#' @param local A logical or list controlling the big data approximation (or \code{NULL} if unresolved)
#'
#' @return A list with elements \code{local} (resolved), \code{obdata},
#'   \code{xcoord}, \code{ycoord}, \code{newdata} (anisotropy-transformed and,
#'   if \code{newdata} is an \code{sf} object, converted to a data frame),
#'   \code{add_newdata_rows}, \code{spcov_params_val}, \code{randcov_params_val},
#'   \code{newdata_model}, and \code{offset} -- everything \code{predict.splm()}
#'   needs after argument matching but before dispatching on \code{interval}/\code{type},
#'   gathered here so \code{predict.splm()} itself does not have to carry this
#'   setup logic inline (its caller stores these directly in its environment
#'   via \code{list2env()}, so none of its downstream code has to change)
#'
#' @noRd
get_prediction_object_splm <- function(object, newdata, scale, local) {
  # check scale is numeric (if specified)
  if (!is.null(scale) && !is.numeric(scale)) {
    stop("scale must be numeric.", call. = FALSE)
  }

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }

  # rename relevant quantities
  obdata <- object$obdata
  xcoord <- object$xcoord
  ycoord <- object$ycoord

  # write newdata if predicting missing data
  if (missing(newdata)) {
    add_newdata_rows <- TRUE
    newdata <- object$newdata
  } else {
    add_newdata_rows <- FALSE
  }

  # edit newdata if random effect or partition factor levels are new
  if (!is.null(object$random) || !is.null(object$partition_factor)) {
    random_names <- all.vars(object$random)
    partition_names <- all.vars(object$partition_factor)
    varnames <- unique(c(random_names, partition_names))
    newdata <- replace_newdata(varnames, obdata, newdata)
  }

  # deal with local
  # exact (non-local) prediction requires inverting an n x n covariance
  # matrix, which becomes impractically slow/memory-hungry past a few
  # thousand observations -- auto-switch to the local approximation above
  # that threshold unless the caller explicitly requested local = FALSE
  if (is.null(local)) {
    if (object$n > 10000) {
      local <- TRUE
      message("Because the sample size of the fitted model object exceeds 10,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun predict() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else {
      local <- FALSE
    }
  }

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save randcov param vector
  randcov_params_val <- coef(object, type = "randcov")

  # sp package classes are unsupported (superseded by sf); detect them via
  # their class's package attribute and point the user to sf instead
  attr_sp <- attr(class(newdata), "package")
  if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
    stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
  }

  if (inherits(newdata, "sf")) {
    # collapse polygons/lines to their centroid coordinates so prediction can
    # proceed on point coordinates like the rest of the package expects
    newdata <- suppressWarnings(sf::st_centroid(newdata))

    newdata <- sf_to_df(newdata)
    names(newdata)[[which(names(newdata) == ".xcoord")]] <- as.character(xcoord) # only relevant if newdata is sf data is not
    names(newdata)[[which(names(newdata) == ".ycoord")]] <- as.character(ycoord) # only relevant if newdata is sf data is not
  }

  # add back in zero column to cover anisotropy (should make anisotropy only available 1-d)
  if (object$dim_coords == 0) {
    # spcov_type "none"/"ie" fit without xcoord/ycoord: both coordinates are
    # synthetic placeholders (see get_point_ref_coords()), so newdata never
    # has them either -- fill both rather than requiring check_newdata_coords()
    # below to find columns that were never meant to exist in newdata
    obdata[[xcoord]] <- 0
    obdata[[ycoord]] <- 0
    newdata[[xcoord]] <- 0
    newdata[[ycoord]] <- 0
  } else if (object$dim_coords == 1) {
    obdata[[ycoord]] <- 0
    newdata[[ycoord]] <- 0
  }

  check_newdata_coords(newdata, xcoord, ycoord)

  if (object$anisotropy) { # could just do rotate != 0 || scale != 1
    # rotate/rescale both observed and new coordinates into the isotropic
    # space the covariance function was fit in, so distances (and therefore
    # predictions) are computed consistently with model fitting
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

  newdata_model_list <- get_newdata_model_matrix(object, newdata)
  newdata <- newdata_model_list$newdata
  newdata_model <- newdata_model_list$newdata_model
  offset <- newdata_model_list$offset
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  # newdata_model may contain columns not present in the fitted model's
  # design matrix (e.g. new/unused factor levels dropped by the fit); keep
  # only the columns that match so newdata_model aligns with betahat
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # finding rows w/out NA
  ob_predictors <- complete.cases(newdata_model)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }

  list(
    local = local, obdata = obdata, xcoord = xcoord, ycoord = ycoord,
    newdata = newdata, add_newdata_rows = add_newdata_rows,
    spcov_params_val = spcov_params_val, randcov_params_val = randcov_params_val,
    newdata_model = newdata_model, offset = offset
  )
}

#' Build the prediction setup object used by \code{predict.spautor()}
#'
#' @param object A fitted model object from [spautor()]
#' @param newdata Data requiring prediction; \code{spautor()} predictions are
#'   always for \code{object}'s missing-index rows, so \code{newdata} is only
#'   used to validate the call (see \code{check_newdata_areal()}) -- either
#'   omitted, or identical to \code{object$newdata}
#' @param scale A scale multiplier for prediction standard errors (or \code{NULL})
#' @param local A logical or list controlling the big data approximation (or \code{NULL} if unresolved)
#'
#' @return A list with elements \code{local} (resolved), \code{newdata},
#'   \code{spcov_params_val}, \code{randcov_params_val}, \code{newdata_model},
#'   and \code{offset} -- everything \code{predict.spautor()} needs after
#'   argument matching but before dispatching on \code{interval}/\code{type},
#'   gathered here so \code{predict.spautor()} itself does not have to carry
#'   this setup logic inline (its caller stores these directly in its
#'   environment via \code{list2env()}, so none of its downstream code has to change)
#'
#' @noRd
get_prediction_object_spautor <- function(object, newdata, scale, local) {
  # deal with local
  if (is.null(local)) {
    local <- FALSE
  }

  # check scale is numeric (if specified)
  if (!is.null(scale) && !is.numeric(scale)) {
    stop("scale must be numeric.", call. = FALSE)
  }

  # error if newdata missing from arguments and object, or if a newdata other
  # than object$newdata was supplied (prediction locations are fixed at
  # fitting time for spautor() -- see check_newdata_areal())
  check_newdata_areal(object$newdata, !missing(newdata), if (missing(newdata)) NULL else newdata, "spautor")

  # write newdata if predicting missing data
  newdata <- object$data[object$missing_index, , drop = FALSE]

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save randcov param vector
  randcov_params_val <- coef(object, type = "randcov")

  newdata_model_list <- get_newdata_model_matrix(object, newdata)
  newdata <- newdata_model_list$newdata
  newdata_model <- newdata_model_list$newdata_model
  offset <- newdata_model_list$offset
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  # newdata_model may contain columns not present in the fitted model's
  # design matrix (e.g. new/unused factor levels dropped by the fit); keep
  # only the columns that match so newdata_model aligns with betahat
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # finding rows w/out NA
  # this isn't really needed, because the error should come on model building
  # but someone could accidentally write over their newdata object after fitting
  # so it is good to keep for completeness
  ob_predictors <- complete.cases(newdata_model)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }

  list(
    local = local, newdata = newdata, spcov_params_val = spcov_params_val,
    randcov_params_val = randcov_params_val, newdata_model = newdata_model, offset = offset
  )
}
