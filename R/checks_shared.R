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
    # should also be given for sf objects
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
