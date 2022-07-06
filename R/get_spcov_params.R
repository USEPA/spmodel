#' Get spatial covariance parameters
#'
#' @param spcov_type The spatial covariance type
#' @param spcov_orig_val Spatial covariance parameters
#'
#' @return A \code{spcov_params} object
#'
#' @noRd
get_spcov_params <- function(spcov_type, spcov_orig_val) {
  if (spcov_type %in% c("exponential", "spherical", "gaussian", "triangular", "circular", "none", "cubic", "penta", "cosine", "wave", "jbessel", "gravity", "rquad", "magnetic")) {
    spcov_params_val <- spcov_params(
      spcov_type = spcov_type,
      de = spcov_orig_val[["de"]],
      ie = spcov_orig_val[["ie"]],
      range = spcov_orig_val[["range"]],
      rotate = spcov_orig_val[["rotate"]],
      scale = spcov_orig_val[["scale"]]
    )
  } else if (spcov_type %in% c("matern", "cauchy", "pexponential")) {
    spcov_params_val <- spcov_params(
      spcov_type = spcov_type,
      de = spcov_orig_val[["de"]],
      ie = spcov_orig_val[["ie"]],
      range = spcov_orig_val[["range"]],
      extra = spcov_orig_val[["extra"]],
      rotate = spcov_orig_val[["rotate"]],
      scale = spcov_orig_val[["scale"]]
    )
  } else if (spcov_type %in% c("car", "sar")) {
    spcov_params_val <- spcov_params(
      spcov_type = spcov_type, de = spcov_orig_val[["de"]],
      ie = spcov_orig_val[["ie"]],
      range = spcov_orig_val[["range"]],
      extra = spcov_orig_val[["extra"]]
    )
  }
  spcov_params_val
}
