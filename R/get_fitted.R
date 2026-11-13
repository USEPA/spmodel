#' Get Fitted Values
#'
#' @param betahat Vector of fixed effects
#' @param spcov_params A \code{spcov_params} object
#' @param X Model matrix
#' @param eigenprods A \code{eigenprods} object
#' @param dist_matrix A distance matrix
#' @param randcov_params A \code{randcov_params} object
#' @param randcov_Zs Random effect design matrices
#'
#' @return A list of fitted values
#'
#' @noRd
get_fitted_splm <- function(betahat, spcov_params, data_object, eigenprods_list,
                            dist_matrix_list, randcov_params = NULL) {
  # X %*% betahat gives the fixed-effects (trend) component of the fit for each
  # partition; these are later combined with the spatial BLUP components below
  fitted_response <- as.numeric(do.call("rbind", lapply(data_object$X_list, function(x) x %*% betahat)))

  if (!is.null(data_object$offset)) {
    fitted_response <- fitted_response + data_object$offset
  }

  ## fitted values
  ### resid pearson is siginv^(-1/2)(y - x beta)
  ### so upchol is from siginv^(-1/2)uptri * siginv^(-1/2)(y - x beta)
  ### packsolve used because upper triangluar
  ### this gives siginv (y - x beta)
  SigInv_r_list <- lapply(eigenprods_list, function(x) x$SigInv_y - x$SigInv_X %*% betahat)

  # covariance matrix built from only the dependent-error (spatially correlated)
  # component -- the independent-error (nugget) component contributes nothing to
  # the BLUP of the spatial signal, since it has zero covariance with any other observation
  spcov_matrix_de_only_list <- lapply(
    dist_matrix_list,
    function(x) spcov_matrix_de(spcov_params = spcov_params, dist_matrix = x)
  )


  ### cov(dependent error) * Z (identity) * siginv (y - x beta)
  fitted_de <- as.numeric(do.call("rbind", mapply(
    s = spcov_matrix_de_only_list, r = SigInv_r_list,
    function(s, r) s %*% r, SIMPLIFY = FALSE
  )))

  ### cov(independent error) is zero so this gives
  ### sigma^2(independent) * Identity * Z (identity) * siginv (y - x beta)
  fitted_ie <- as.numeric(spcov_params[["ie"]] * do.call("rbind", SigInv_r_list))

  ## fitted random effects
  fitted_randcov <- get_fitted_randcov(data_object, randcov_params, SigInv_r_list)

  fitted_values <- list(
    response = fitted_response,
    spcov = list(de = fitted_de, ie = fitted_ie),
    randcov = fitted_randcov
  )
}

#' Get fitted values for an \code{spautor()} model
#'
#' @param betahat Vector of fixed effects
#' @param spcov_params A \code{spcov_params} object
#' @param data_object The data object
#' @param eigenprods A \code{eigenprods} object
#' @param randcov_params A \code{randcov_params} object
#'
#' @return A list of fitted values
#'
#' @noRd
get_fitted_spautor <- function(betahat, spcov_params, data_object, eigenprods,
                               randcov_params = NULL) {
  # subset the full neighborhood/weights structures down to just the observed
  # locations (spautor models can include unobserved locations for prediction)
  dist_matrix <- data_object$W[data_object$observed_index, data_object$observed_index, drop = FALSE]
  M <- data_object$M[data_object$observed_index]

  fitted_response <- data_object$X %*% betahat

  if (!is.null(data_object$offset)) {
    fitted_response <- fitted_response + data_object$offset
  }

  ## fitted values
  ### resid pearson is siginv^(-1/2)(y - x beta)
  ### so upchol is from siginv^(-1/2)uptri * siginv^(-1/2)(y - x beta)
  ### packsolve used because upper triangluar
  ### this gives siginv (y - x beta)
  SigInv_r <- eigenprods$SigInv_y - eigenprods$SigInv_X %*% betahat

  # Build the dependent-error covariance without adding any nugget floor.
  spcov_matrix_de_only <- spcov_matrix_de(spcov_params = spcov_params, dist_matrix = dist_matrix, M = M)

  if (!is.null(data_object$partition_factor)) {
    # zero out covariance between observations in different partitions -- a
    # partition factor means observations are only allowed to be
    # correlated with others sharing the same partition level
    spcov_matrix_de_only <- spcov_matrix_de_only * data_object$partition_matrix[data_object$observed_index, data_object$observed_index, drop = FALSE]
  }

  ### cov(dependent error) * Z (identity) * siginv (y - x beta)
  fitted_de <- spcov_matrix_de_only %*% SigInv_r

  ### cov(independent error) is zero so this gives
  ### sigma^2(independent) * Identity * Z (identity) * siginv (y - x beta)
  fitted_ie <- spcov_params[["ie"]] * SigInv_r

  ## fitted random effects
  fitted_randcov <- get_fitted_randcov(data_object, randcov_params, list(SigInv_r))

  fitted_values <- list(
    response = as.numeric(fitted_response),
    spcov = list(de = as.numeric(fitted_de), ie = as.numeric(fitted_ie)),
    randcov = fitted_randcov
  )
}
