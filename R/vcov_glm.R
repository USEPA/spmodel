#' @rdname vcov.spmodel
#' @method vcov spglm
#' @order 3
#' @export
vcov.spglm <- function(object, var_correct = TRUE, ...) {
  # type is hard-coded (see vcov.splm) since only fixed effects are supported
  type <- "fixed"
  if (type == "fixed") {
    # "corrected" adjusts the naive fixed-effect covariance to account for the
    # extra uncertainty from estimating the latent random effects (see
    # get_wts_varw()); "uncorrected" is that naive (asymptotic) estimate
    if (var_correct) {
      return(object$vcov$fixed$corrected)
    } else {
      return(object$vcov$fixed$uncorrected)
    }
  }
}

#' @rdname vcov.spmodel
#' @method vcov spgautor
#' @order 4
#' @export
vcov.spgautor <- vcov.spglm
