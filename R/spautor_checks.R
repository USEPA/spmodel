#' Perform initial checks for spautor models
#'
#' @param spcov_type The spatial covariance type
#' @param W_given Is the spatial weight matrix given?
#' @param data data
#'
#' @return Error messages or nothing
#'
#' @noRd
spautor_checks <- function(spcov_type, W_given, data, estmethod) {

  # if (!requireNamespace("sf", quietly = TRUE) && !W_given) {
  #   stop("When W is not provided, package \"sf\" must be installed to use spautor().",
  #     call. = FALSE
  #   )
  # }

  if (spcov_type %in% c(
    "exponential", "spherical", "gaussian", "triangular",
    "circular", "cubic", "pentaspherical", "cosine", "wave",
    "jbessel", "gravity", "rquad", "magnetic",
    "matern", "cauchy", "pexponential", "none"
  )) {
    stop("Invalid spatial covariance type for spautor(). To fit models for point-referenced data, use splm().", call. = FALSE)
  }

  # return an error if data are not the correct spcov_type
  if (!W_given && !inherits(data, c("SpatialPolygonsDataFrame", "sf"))) {
    stop("Data must be a SpatialPolygonsDataFrame (sp object) or an sf object", call. = FALSE)
  }

  if (!estmethod %in% c("reml", "ml")) {
    stop("Estimation method must be \"reml\", or \"ml\".", call. = FALSE)
  }
}
