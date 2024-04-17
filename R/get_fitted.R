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

  # cov params no de   (set ie portion to zero because BLUP only uses cov(dependent error))
  spcov_params_de_only <- spcov_params
  spcov_params_de_only[["ie"]] <- 0
  spcov_matrix_de_only_list <- lapply(
    dist_matrix_list,
    function(x) spcov_matrix(spcov_params = spcov_params_de_only, dist_matrix = x)
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
  if (is.null(names(randcov_params))) {
    fitted_randcov <- NULL
  } else {
    fitted_randcov <- lapply(names(randcov_params), function(x) {
      fitted_val <- randcov_params[[x]] * do.call("rbind", mapply(
        z = data_object$randcov_list,
        r = SigInv_r_list,
        function(z, r) {
          crossprod(z[[x]][["Z"]], r)
        }
      ))
      fitted_val <- tapply(fitted_val, rownames(fitted_val), function(x) {

        if (any(x != 0)) {
          val <- mean(x[x != 0])
        } else {
          val <- 0
        }
        # if (length(val) == 0) { # replace if all zeros somehow
        #   val <- rep(0, length(x))
        #   names(val) <- names(x)
        # }
        # if all elements of x are zero then val is NaN, so reset it to zero
        # this works but above uncommented code is cleaner
        # if (is.na(val)) {
        #   val <- 0
        # }
        # val
      })
      # all combinations yields values with many zeros -- don't want to include these in the mean
      names_fitted_val <- rownames(fitted_val)
      fitted_val <- as.numeric(fitted_val)
      names(fitted_val) <- names_fitted_val
      fitted_val
    })
    names(fitted_randcov) <- names(randcov_params)
  }

  fitted_values <- list(
    response = fitted_response,
    spcov = list(de = fitted_de, ie = fitted_ie),
    randcov = fitted_randcov
  )
}

get_fitted_spautor <- function(betahat, spcov_params, data_object, eigenprods,
                               randcov_params = NULL) {
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

  # cov params no de   (set ie portion to zero because BLUP only uses cov(dependent error))
  spcov_params_de_only <- spcov_params
  spcov_params_de_only[["ie"]] <- 0
  spcov_matrix_de_only <- spcov_matrix(spcov_params = spcov_params_de_only, dist_matrix = dist_matrix, M = M)

  if (!is.null(data_object$partition_factor)) {
    spcov_matrix_de_only <- spcov_matrix_de_only * data_object$partition_matrix[data_object$observed_index, data_object$observed_index, drop = FALSE]
  }

  ### cov(dependent error) * Z (identity) * siginv (y - x beta)
  fitted_de <- spcov_matrix_de_only %*% SigInv_r

  ### cov(independent error) is zero so this gives
  ### sigma^2(independent) * Identity * Z (identity) * siginv (y - x beta)
  fitted_ie <- spcov_params[["ie"]] * SigInv_r

  ## fitted random effects
  if (is.null(names(randcov_params))) {
    fitted_randcov <- NULL
  } else {
    if (is.null(data_object$partition_factor)) {
      ob_randcov_Zs <- get_randcov_Zs(data_object$obdata, names(randcov_params), ZZt = FALSE)
      fitted_randcov <- lapply(names(randcov_params), function(x) {
        fitted_val <- randcov_params[[x]] * crossprod(ob_randcov_Zs[[x]][["Z"]], SigInv_r)
        names_fitted_val <- rownames(fitted_val)
        fitted_val <- as.vector(fitted_val)
        names(fitted_val) <- names_fitted_val
        fitted_val
      })
      names(fitted_randcov) <- names(randcov_params)
    } else {
      index <- unname(model.response(model.frame(reformulate("1", response = labels(terms(data_object$partition_factor))),
        data = data_object$obdata
      )))
      index_val <- unique(index)
      ob_randcov_Zs <- get_randcov_Zs(data_object$obdata, names(randcov_params), ZZt = FALSE)
      fitted_randcov <- lapply(names(randcov_params), function(x) {
        fitted_val <- lapply(index_val, function(y) {
          row_val <- y == index
          fitted_vals <- randcov_params[[x]] *
            crossprod(ob_randcov_Zs[[x]][["Z"]][row_val, , drop = FALSE], SigInv_r[row_val, , drop = FALSE])
        })
        fitted_val <- do.call("rbind", fitted_val)
        fitted_val <- tapply(fitted_val, rownames(fitted_val), function(x) {
          val <- mean(x[x != 0])
          if (length(val) == 0) { # replace if all zeros somehow
            val <- rep(0, length(x))
            names(val) <- names(x)
          }
          val
        })
        names_fitted_val <- rownames(fitted_val)
        fitted_val <- as.vector(fitted_val)
        names(fitted_val) <- names_fitted_val
        fitted_val
      })
      names(fitted_randcov) <- names(randcov_params)
    }
  }




  fitted_values <- list(
    response = as.numeric(fitted_response),
    spcov = list(de = as.numeric(fitted_de), ie = as.numeric(fitted_ie)),
    randcov = fitted_randcov
  )
}
