#' Build the prediction setup object used by \code{predict.spglm()}
#'
#' @param object A fitted model object from [spglm()]
#' @param newdata Data requiring prediction (possibly missing, to predict \code{object$newdata})
#' @param dispersion A dispersion value overriding \code{object}'s fitted
#'   dispersion (or \code{NULL} to use the fitted value)
#' @param newdata_size Binomial trial sizes for \code{newdata} (resolved to
#'   \code{NULL} by the caller if not supplied)
#' @param local A logical or list controlling the big data approximation (or \code{NULL} if unresolved)
#'
#' @return A list with elements \code{object} (with \code{dispersion} applied,
#'   if supplied), \code{local} (resolved), \code{newdata_size} (resolved),
#'   \code{obdata}, \code{xcoord}, \code{ycoord}, \code{newdata}, \code{add_newdata_rows},
#'   \code{spcov_params_val}, \code{dispersion_params_val}, \code{randcov_params_val},
#'   \code{newdata_model}, and \code{offset} -- everything \code{predict.spglm()}
#'   needs after argument matching but before dispatching on \code{interval}/\code{type},
#'   gathered here so \code{predict.spglm()} itself does not have to carry this
#'   setup logic inline (its caller stores these directly in its environment
#'   via \code{list2env()}, so none of its downstream code has to change)
#'
#' @noRd
get_prediction_object_spglm <- function(object, newdata, dispersion, newdata_size, local) {
  # handle dispersion argument if provided
  # binomial/poisson have no free dispersion parameter (it's fixed at 1 by
  # definition of those families), so overriding it is only meaningful, and
  # only allowed, for families that do have one (e.g. Gaussian, negative binomial)
  if (!is.null(dispersion)) {
    if (object$family %in% c("binomial", "poisson") && dispersion != 1) {
      stop("dispersion is fixed at one for binomial and poisson families.", call. = FALSE)
    }
    object$coefficients$dispersion[1] <- dispersion
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

  if (!is.null(object$random) || !is.null(object$partition_factor)) {
    random_names <- all.vars(object$random)
    partition_names <- all.vars(object$partition_factor)
    varnames <- unique(c(random_names, partition_names))
    newdata <- replace_newdata(varnames, obdata, newdata)
  }

  # set newdata_size if needed
  # binomial predictions need a number-of-trials denominator; if the caller
  # didn't supply one, default to Bernoulli trials (size 1 per row)
  if (is.null(newdata_size) && object$family == "binomial") {
    newdata_size <- rep(1, NROW(newdata))
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

  # save dispersion param vector
  dispersion_params_val <- as.vector(coef(object, type = "dispersion")) # remove class

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
  if (object$dim_coords == 1) {
    obdata[[ycoord]] <- 0
    newdata[[ycoord]] <- 0
  }

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

  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  # workaround: duplicate the single row so poly() has enough rows to build
  # its basis without erroring, build the model matrix, then keep only row 1
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
    # find offset
    offset <- model.offset(newdata_model_frame)
  }
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
    object = object, local = local, newdata_size = newdata_size,
    obdata = obdata, xcoord = xcoord, ycoord = ycoord,
    newdata = newdata, add_newdata_rows = add_newdata_rows,
    spcov_params_val = spcov_params_val, dispersion_params_val = dispersion_params_val,
    randcov_params_val = randcov_params_val, newdata_model = newdata_model, offset = offset
  )
}

#' Build the prediction setup object used by \code{predict.spgautor()}
#'
#' @param object A fitted model object from [spgautor()]
#' @param newdata Data requiring prediction (unused directly; \code{spgautor()}
#'   predictions are always for \code{object}'s missing-index rows -- present
#'   only so \code{missing(newdata)} can still be checked against \code{object$newdata})
#' @param dispersion A dispersion value overriding \code{object}'s fitted
#'   dispersion (or \code{NULL} to use the fitted value)
#' @param newdata_size Binomial trial sizes for \code{newdata} (resolved to
#'   \code{NULL} by the caller if not supplied)
#' @param local A logical or list controlling the big data approximation (or \code{NULL} if unresolved)
#'
#' @return A list with elements \code{object} (with \code{dispersion} applied,
#'   if supplied), \code{local} (resolved), \code{newdata_size} (resolved),
#'   \code{newdata}, \code{spcov_params_val}, \code{dispersion_params_val},
#'   \code{randcov_params_val}, \code{newdata_model}, and \code{offset} --
#'   everything \code{predict.spgautor()} needs after argument matching but
#'   before dispatching on \code{interval}/\code{type}, gathered here so
#'   \code{predict.spgautor()} itself does not have to carry this setup logic
#'   inline (its caller stores these directly in its environment via
#'   \code{list2env()}, so none of its downstream code has to change)
#'
#' @noRd
get_prediction_object_spgautor <- function(object, newdata, dispersion, newdata_size, local) {
  # handle dispersion argument if provided
  # binomial/poisson have no free dispersion parameter (it's fixed at 1 by
  # definition of those families), so overriding it is only meaningful, and
  # only allowed, for families that do have one (e.g. Gaussian, negative binomial)
  if (!is.null(dispersion)) {
    if (object$family %in% c("binomial", "poisson") && dispersion != 1) {
      stop("dispersion is fixed at one for binomial and poisson families.", call. = FALSE)
    }
    object$coefficients$dispersion[1] <- dispersion
  }

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }

  # deal with local
  if (is.null(local)) {
    local <- FALSE
  }

  # write newdata if predicting missing data
  newdata <- object$data[object$missing_index, , drop = FALSE]

  # set newdata_size if needed
  # binomial predictions need a number-of-trials denominator; if the caller
  # didn't supply one, default to Bernoulli trials (size 1 per row)
  if (is.null(newdata_size) && object$family == "binomial") {
    newdata_size <- rep(1, NROW(newdata))
  }

  # save spcov param vector
  spcov_params_val <- coef(object, type = "spcov")

  # save dispersion param vector
  dispersion_params_val <- as.vector(coef(object, type = "dispersion")) # remove class

  # save randcov param vector
  randcov_params_val <- coef(object, type = "randcov")

  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  # workaround: duplicate the single row so poly() has enough rows to build
  # its basis without erroring, build the model matrix, then keep only row 1
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
    # find offset
    offset <- model.offset(newdata_model_frame)
  }
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
    object = object, local = local, newdata_size = newdata_size,
    newdata = newdata, spcov_params_val = spcov_params_val,
    dispersion_params_val = dispersion_params_val, randcov_params_val = randcov_params_val,
    newdata_model = newdata_model, offset = offset
  )
}
