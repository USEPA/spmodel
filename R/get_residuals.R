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
  residuals_response <- as.numeric(do.call("rbind", mapply(
    y = data_object$y_list, x = data_object$X_list,
    function(y, x) y - x %*% betahat, SIMPLIFY = FALSE
  )))
  residuals_pearson <- as.numeric(do.call(
    "rbind",
    lapply(eigenprods_list, function(x) x$SqrtSigInv_y - x$SqrtSigInv_X %*% betahat)
  ))
  residuals_standardized <- residuals_pearson / sqrt(1 - hatvalues) # (I - H on bottom)
  list(response = residuals_response, pearson = residuals_pearson, standardized = residuals_standardized)
}

get_residuals_spautor <- function(betahat, X, y, eigenprods, hatvalues) {
  residuals_response <- as.numeric(y - X %*% betahat)
  residuals_pearson <- as.numeric(eigenprods$SqrtSigInv_y - eigenprods$SqrtSigInv_X %*% betahat)
  residuals_standardized <- residuals_pearson / sqrt(1 - hatvalues) # (I - H on bottom)
  list(response = residuals_response, pearson = residuals_pearson, standardized = residuals_standardized)
}
