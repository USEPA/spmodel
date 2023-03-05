#' Perform initial checks for spautor models
#'
#' @param spcov_type The spatial covariance type
#' @param y_coord_given Is the y-coordinate given?
#'
#' @return Error messages or nothing
#'
#' @noRd
spglm_checks <- function(spcov_initial, xcoord_given, ycoord_given, estmethod, anisotropy, random_given) {
  spcov_type <- class(spcov_initial)
  if (spcov_type %in% c("car", "sar")) {
    stop("Invalid spatial covariance type for spglm(). To fit models for autoregressive data, use spgautor().", call. = FALSE)
  }

  if (spcov_type %in% c("triangular", "cosine") && ycoord_given) {
    warning(paste0(spcov_type, " covariance can only be used in one dimension. Ignoring y-coordinate."), call. = FALSE)
    # should also be given for sf objects
  }

  if (!estmethod %in% c("reml", "ml")) {
    stop("Estimation method must be \"reml\" or \"ml\".", call. = FALSE)
  }
}
