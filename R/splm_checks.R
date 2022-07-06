#' Perform initial checks for spautor models
#'
#' @param spcov_type The spatial covariance type
#' @param y_coord_given Is the y-coordinate given?
#'
#' @return Error messages or nothing
#'
#' @noRd
splm_checks <- function(spcov_initial, xcoord_given, ycoord_given, estmethod, anisotropy, random_given) {
  spcov_type <- class(spcov_initial)
  if (spcov_type %in% c("car", "sar")) {
    stop("Invalid spatial covariance type for splm(). To fit models for autoregressive data, use spautor().", call. = FALSE)
  }
  # } else if (!spcov_type %in% c("exponential", "spherical", "gaussian", "triangular",
  #                               "circular", "cubic", "penta", "cosine", "wave",
  #                               "jbessel", "gravity", "rquad", "magnetic",
  #                               "matern", "cauchy", "pexponential", "none")) {
  #   stop("Invalid spatial covariance type. Valid spatial covariance types include \"exponential\", \"spherical\", \"gaussian\", \"triangular\", \"circular\", \"cubic\", \"penta\", \"cosine\", \"wave\", \"jbessel\", \"gravity\", \"rquad\", \"magnetic\", \"matern\", \"cauchy\", \"pexponential\", and \"none\".")
  # }

  if (spcov_type %in% c("triangular", "cosine") && ycoord_given) {
    warning(paste0(spcov_type, " covariance can only be used in one dimension. Overwriting y-coordinate."), call. = FALSE)
    # should also be given for sf objects
  }

  # if (!xcoord_given && spcov_type != "none" && ) {
  #   stop("The xcoord argument must be specified.", call. = FALSE)
  # }

  if (!estmethod %in% c("reml", "ml", "sv-wls", "sv-cl")) {
    stop("Estimation method must be \"reml\", \"ml\", \"sv-wls\", or \"sv-cl\".", call. = FALSE)
  }

  if (estmethod %in% c("sv-wls", "sv-cl")) {
    if (anisotropy) {
      stop("Anisotropy cannot be estimated if estmethod is \"sv-wls\" or \"sv-cl\". To incorporate known anisotropy parameters, use the spcov_initial argument.", call. = FALSE)
    }
    if (("rotate" %in% names(spcov_initial$is_known) && !spcov_initial$is_known["rotate"]) || ("scale" %in% names(spcov_initial$is_known) && !spcov_initial$is_known["scale"])) {
      stop("When estmethod is  \"sv-wls\" or \"sv-cl\", the anisotropy parameters rotate and scale must be assumed known.")
    }
    if (random_given) {
      stop("Random effects cannot be estimated when estmethod is \"sv-wls\" or \"sv-cl\".", call. = FALSE)
    }
  }
}
