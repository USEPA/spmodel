#' Reject areal (car/sar) covariance types for point-referenced model constructors
#'
#' @param spcov_type The spatial covariance type
#' @param current_fun Name of the calling constructor (e.g. "splm"), used in the error message
#' @param target_fun Name of the areal-data constructor to redirect the user to (e.g. "spautor")
#'
#' @return Error message or nothing
#'
#' @noRd
check_not_areal_type <- function(spcov_type, current_fun, target_fun) {
  # car/sar are areal (autoregressive) covariance types and belong to the areal
  # constructor, not the point-referenced one calling this check
  if (spcov_type %in% c("car", "sar")) {
    stop(
      paste0(
        "Invalid spatial covariance type for ", current_fun, "(). To fit models for autoregressive data, use ",
        target_fun, "()."
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Warn that a one-dimensional-only covariance type ignores a supplied y-coordinate
#'
#' @param spcov_type The spatial covariance type
#' @param ycoord_given Is the y-coordinate given?
#'
#' @return Warning or nothing
#'
#' @noRd
warn_ycoord_ignored_for_1d_cov <- function(spcov_type, ycoord_given) {
  if (spcov_type %in% c("triangular", "cosine") && ycoord_given) {
    warning(paste0(spcov_type, " covariance can only be used in one dimension. Ignoring y-coordinate."), call. = FALSE)
  }
  invisible(NULL)
}

#' Reject point-referenced (geostatistical) covariance types for areal model constructors
#'
#' @param spcov_type The spatial covariance type
#'
#' @return Error message or nothing
#'
#' @noRd
check_not_point_referenced_type <- function(spcov_type) {
  # these covariance types are for point-referenced (geostatistical) data and
  # are handled by splm()/spglm(); the areal constructors only support car/sar
  if (spcov_type %in% c(
    "exponential", "spherical", "gaussian", "triangular",
    "circular", "cubic", "pentaspherical", "cosine", "wave",
    "jbessel", "gravity", "rquad", "magnetic",
    "matern", "cauchy", "pexponential", "none", "ie"
  )) {
    stop("Invalid spatial covariance type for spautor(). To fit models for point-referenced data, use splm().", call. = FALSE)
  }
  invisible(NULL)
}

#' Require a spatial weight matrix or a data object that can build one
#'
#' @param W_given Is the spatial weight matrix given?
#' @param data data
#'
#' @return Error message or nothing
#'
#' @noRd
check_W_given_data_class <- function(W_given, data) {
  # return an error if data are not the correct spcov_type
  if (!W_given && !inherits(data, c("SpatialPolygonsDataFrame", "sf"))) {
    stop("Data must be a SpatialPolygonsDataFrame (sp object) or an sf object", call. = FALSE)
  }
  invisible(NULL)
}

#' Require estmethod to be "reml" or "ml"
#'
#' @param estmethod The estimation method
#' @param message_text The exact error message to use (callers' wording differs slightly)
#'
#' @return Error message or nothing
#'
#' @noRd
check_estmethod_reml_ml <- function(estmethod, message_text) {
  if (!estmethod %in% c("reml", "ml")) {
    stop(message_text, call. = FALSE)
  }
  invisible(NULL)
}

#' Require family to be one of spmodel's supported GLM families
#'
#' @param family The response family
#'
#' @return Error message or nothing
#'
#' @noRd
check_family_valid <- function(family) {
  family_valid <- c("binomial", "poisson", "nbinomial", "Gamma", "inverse.gaussian", "beta")
  if (!(family %in% family_valid)) {
    stop(paste(family, " is not a valid glm family.", sep = ""), call. = FALSE)
  }
  invisible(NULL)
}

#' Require every variable used in formula/random/partition_factor to exist in data
#'
#' @param formula The model formula
#' @param data The data
#' @param random An optional random effect formula (or \code{NULL})
#' @param partition_factor An optional partition factor formula (or \code{NULL})
#'
#' @details Without this check, a variable used in \code{formula}/\code{random}/
#'   \code{partition_factor} but absent from \code{data} is not caught here.
#'
#' @return Error message or nothing
#'
#' @noRd
check_formula_vars_in_data <- function(formula, data, random = NULL, partition_factor = NULL) {
  if ("." %in% all.vars(random)) {
    stop("The `.` shorthand is not supported in random. Explicitly list the desired variable(s).", call. = FALSE)
  }
  if ("." %in% all.vars(partition_factor)) {
    stop("The `.` shorthand is not supported in partition_factor. Explicitly list the desired variable(s).", call. = FALSE)
  }
  formula_vars <- unique(c(all.vars(formula), all.vars(random), all.vars(partition_factor)))
  formula_vars <- setdiff(formula_vars, ".")
  missing_vars <- setdiff(formula_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "Variable(s) ", paste0("\"", missing_vars, "\"", collapse = ", "),
      " used in formula, random, or partition_factor not found in data.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Reject newdata for areal (spautor()/spgautor()) prediction unless it's object$newdata
#'
#' @param object_newdata The fitted model's own \code{$newdata} (or \code{NULL}
#'   if there were no missing-response rows at fitting time)
#' @param newdata_given Whether the caller's \code{newdata} argument to
#'   \code{predict()} was supplied (i.e. \code{!missing(newdata)})
#' @param newdata The caller's \code{newdata} value; only accessed when
#'   \code{newdata_given} is \code{TRUE} (evaluating a genuinely missing
#'   argument would itself error)
#' @param current_fun Name of the calling constructor, used in the error
#'   message (e.g. \code{"spautor"})
#'
#' @details Unlike \code{splm()}/\code{spglm()}, \code{spautor()}/\code{spgautor()}
#'   prediction locations are fixed when the model is fit, as they determine
#'   the neighbor structure (\code{W}/\code{M}) used throughout fitting. This
#'   implies \code{newdata} cannot be provided by the user if it is different
#'   from object$newdata.  
#'
#' @return Error message or nothing
#'
#' @noRd
check_newdata_areal <- function(object_newdata, newdata_given, newdata, current_fun) {
  if (!newdata_given) {
    if (is.null(object_newdata)) {
      stop("No missing data to predict. Fit the model with NA response values for the locations you want to predict.", call. = FALSE)
    }
  } else if (!identical(newdata, object_newdata)) {
    stop(
      "newdata cannot be specified for ", current_fun, "() model objects different from object$newdata, ",
      "because prediction locations are fixed when the model is fit (they determine the neighbor structure used in fitting). ",
      "Ignoring newdata and predicting for object$newdata instead.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Resolve \code{newdata} for \code{augment()} on an autoregressive model
#'
#' @param object A fitted \code{spautor()}/\code{spgautor()} model object
#' @param se_fit The \code{se_fit} argument to \code{augment}
#' @param interval The \code{interval} argument to \code{augment}
#' 
#' @details \code{se_fit} and \code{interval} describe predictions, and for
#'   \code{spautor()}/\code{spgautor()} predictions exist only at the locations
#'   set aside when the model was fit. Those locations are part of the model:
#'   they enter the neighbor structure (\code{W}/\code{M}) and so determine the
#'   covariance of the observed data itself, which means a standard error at an
#'   observed location would imply a (new) observation that is different from the model than the one that was
#'   fit. \code{object$newdata} is therefore the only value \code{newdata} can
#'   take, so \code{augment()}
#'   supplies it leaves a message, because the returned rows change from the
#'   observed data to the prediction locations.
#'
#'   \code{splm()}/\code{spglm()} have no such restriction: their covariance
#'   does not depend on where predictions are made, so \code{se_fit} and
#'   \code{interval = "confidence"} are well defined at the observed locations
#'   and are computed there directly.
#'
#' @return \code{object$newdata} when \code{se_fit} or \code{interval} asks for
#'   a prediction quantity, otherwise \code{NULL}
#'
#' @noRd
augment_areal_newdata <- function(object, se_fit, interval) {
  if (!se_fit && interval == "none") {
    return(NULL)
  }
  if (is.null(object$newdata)) {
    stop(
      "No missing data to predict. Fit the model with NA response values for the locations you want to predict.",
      call. = FALSE
    )
  }
  message(
    "se_fit and interval describe predictions, which for autoregressive models are only defined at the ",
    "locations set aside when the model was fit (they determine the neighbor structure, and hence the ",
    "covariance of the observed data). Returning augmented output for object$newdata rather than the observed data."
  )
  object$newdata
}

#' Warn when \code{augment()} is given an \code{interval} it cannot use
#'
#' @param interval The caller's (already matched) \code{interval} argument
#' @param newdata_given Whether \code{newdata} was supplied to \code{augment()}
#'
#' @details Without \code{newdata}, \code{augment()} describes the observed
#'   data. A confidence interval is still meaningful there -- it brackets the
#'   fitted mean -- but a prediction interval is not, because it needs a
#'   location that has not been observed. \code{interval = "prediction"} is
#'   therefore downgraded to \code{"none"} with a warning rather than silently
#'   producing a degenerate interval.
#'
#' @return The \code{interval} to use, downgraded to \code{"none"} if it cannot
#'   be honored
#'
#' @noRd
check_interval_augment <- function(interval, newdata_given) {
  if (!newdata_given && interval == "prediction") {
    warning(
      "interval = \"prediction\" is ignored when newdata is not supplied, because a prediction interval ",
      "requires a location that has not been observed. Supply newdata for prediction intervals, or use ",
      "interval = \"confidence\" for an interval around the fitted mean at the observed locations.",
      call. = FALSE
    )
    interval <- "none"
  }
  interval
}

