#' Fill in default (\code{NA}, unknown) spatial covariance initial values for GLM-type models
#'
#' @param family The response family
#' @param spcov_initial A \code{spcov_initial} object
#' @param anisotropy Whether anisotropy is being modeled
#' @param is_W_connected For areal models, whether the neighbor graph is fully connected
#'
#' @return A \code{spcov_initial} object with defaults filled in, as
#'   \code{spcov_initial_NA()} does for \code{splm()}/\code{spautor()}, plus
#'   fixing \code{ie} at zero for the \code{"none"} covariance type, since
#'   GLM-type models without spatial dependence have no need for a separate
#'   independent error term beyond the dispersion parameter
#'
#' @noRd
spcov_initial_NA_glm <- function(family, spcov_initial, anisotropy = FALSE, is_W_connected = NULL) {
  spcov_initial_NA_val <- spcov_initial_NA(spcov_initial, anisotropy, is_W_connected)
  spcov_initial_NA_glm_val <- spcov_initial_NA_val

  # unlike splm()/spautor(), GLM-type models with no spatial covariance don't
  # need an estimable independent error variance (dispersion already covers
  # that role), so ie is overridden to a fixed, known 0 here
  if (inherits(spcov_initial_NA_glm_val, "none")) {
    spcov_initial_NA_glm_val$initial["ie"] <- 0
    spcov_initial_NA_glm_val$is_known["ie"] <- TRUE
  }

  spcov_initial_NA_glm_val
}
