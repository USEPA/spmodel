#' Get residuals
#'
#' @param betahat Fixed effects
#' @param X Model matrix
#' @param y Response vector
#' @param eigenprods A \code{eigenprods} object
#' @param hatvalues Leverage values
#'
#' @return A list of relevant residuals
#'
#' @noRd
get_residuals_splm <- function(betahat, data_object, eigenprods_list, hatvalues) {
  # data is stored split by group (e.g. partition factor levels), so residuals
  # are computed per group with mapply() and then stacked back into one vector
  residuals_response <- as.numeric(do.call("rbind", mapply(
    y = data_object$y_list, x = data_object$X_list,
    function(y, x) y - x %*% betahat, SIMPLIFY = FALSE
  )))
  # Pearson (whitened) residuals: y and X are already pre-multiplied by the
  # inverse-covariance square root (SqrtSigInv), so this is (Sigma^-1/2)(y - X betahat)
  residuals_pearson <- as.numeric(do.call(
    "rbind",
    lapply(eigenprods_list, function(x) x$SqrtSigInv_y - x$SqrtSigInv_X %*% betahat)
  ))
  residuals_standardized <- residuals_pearson / sqrt(1 - hatvalues) # (I - H on bottom)
  list(response = residuals_response, pearson = residuals_pearson, standardized = residuals_standardized)
}

#' Get residuals for an \code{spautor()} model
#'
#' @param betahat Fixed effects
#' @param X Model matrix
#' @param y Response vector
#' @param eigenprods A \code{eigenprods} object
#' @param hatvalues Leverage values
#'
#' @return A list of relevant residuals
#'
#' @noRd
get_residuals_spautor <- function(betahat, X, y, eigenprods, hatvalues) {
  residuals_response <- as.numeric(y - X %*% betahat)
  # Pearson (whitened) residuals using the pre-whitened y/X (Sigma^-1/2 applied already)
  residuals_pearson <- as.numeric(eigenprods$SqrtSigInv_y - eigenprods$SqrtSigInv_X %*% betahat)
  residuals_standardized <- residuals_pearson / sqrt(1 - hatvalues) # (I - H on bottom)
  list(response = residuals_response, pearson = residuals_pearson, standardized = residuals_standardized)
}
