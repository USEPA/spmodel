#' Adjust the covariance matrix of the fixed effects
#'
#' @param invcov_betahat_list A list of inverse covariances of beta hat for each grouping
#' @param betahat_list A list of beta hat for each grouping
#' @param betahat Average betahat
#' @param cholprods_list A list of cholesky products for each grouping
#' @param data_object The data object
#' @param spcov_params Spatial covariance parameters
#' @param randcov_params Random effects
#' @param cov_betahat_noadjust The covariance of betahat with no adjustment
#' @param var_adjust The variance adjustment type. Options include \code{"none"},
#'   \code{"theoretical"}, \code{"empirical"}, and \code{"pooled"}.
#'
#' @return An adjusted covariance matrix of betahat
#' @noRd
cov_betahat_adjust <- function(invcov_betahat_list, betahat_list,
                               betahat, cholprods_list, data_object, spcov_params,
                               randcov_params, cov_betahat_noadjust, var_adjust) {
  P <- length(betahat_list)
  # reset var_adjust if only one partition
  if (P == 1 || inherits(spcov_params, "none")) {
    var_adjust <- "none"
  }
  # reset var_adjust if partitioning used but no local option used
  if (!is.null(data_object$partition_factor_initial) && is.null(data_object$partition_factor)) {
    var_adjust <- "none"
  }

  if (var_adjust == "empirical") {
    cov_betahat_adjust_list <- lapply(betahat_list, function(x) {
      betahat_diff <- x - betahat
      tcrossprod(betahat_diff, betahat_diff)
    })
    cov_betahat_adjust_val <- Reduce("+", cov_betahat_adjust_list) / (P * (P - 1))
  } else if (var_adjust == "pooled") {
    cov_betahat_list <- tryCatch(
      error = function(cnd) NULL,
      lapply(invcov_betahat_list, function(x) chol2inv(chol(forceSymmetric(x))))
    )
    if (is.null(cov_betahat_list)) {
      var_adjust <- "none"
      warning(
        "At least one partition's inverse covariance matrix is singular. Redjusting using var_adjust = \"none\".",
        call. = FALSE
      )
    } else {
      cov_betahat_adjust_val <- Reduce("+", cov_betahat_list) / P^2
    }
  }

  if (var_adjust == "theoretical") {
    index_grid <- expand.grid(d1 = seq_len(P), d2 = seq_len(P))
    index_grid <- index_grid[index_grid$d2 > index_grid$d1, , drop = FALSE]
    index_list <- split(index_grid, seq_len(NROW(index_grid)))
    if (data_object$parallel) {
      W_adjust_list <- parallel::parLapply(data_object$cl, index_list, get_W_ij_parallel,
        cholprods_list = cholprods_list,
        spcov_params = spcov_params, randcov_params = randcov_params,
        randcov_names = names(randcov_params), data_object = data_object
      )
    } else {
      W_adjust_list <- lapply(index_list, function(x) {
        get_W_ij(
          x$d1, x$d2, cholprods_list,
          spcov_params, randcov_params,
          names(randcov_params), data_object
        )
      })
    }
    W_adjust <- Reduce("+", W_adjust_list)
    cov_betahat_adjust_val <- cov_betahat_noadjust + cov_betahat_noadjust %*% W_adjust %*% cov_betahat_noadjust
  } else if (var_adjust == "none") {
    cov_betahat_adjust_val <- cov_betahat_noadjust
  }

  cov_betahat_adjust_val
}

get_W_ij <- function(d1_index, d2_index, cholprods_list, spcov_params,
                     randcov_params, randcov_names, data_object) {
  d1 <- data_object$obdata_list[[d1_index]]
  d2 <- data_object$obdata_list[[d2_index]]

  # anisotropy correction
  if (data_object$anisotropy) {
    anis_d1 <- transform_anis(
      data = d1,
      data_object$xcoord, data_object$ycoord,
      spcov_params[["rotate"]], spcov_params[["scale"]]
    )
    d1[[data_object$xcoord]] <- anis_d1$xcoord_val
    d1[[data_object$ycoord]] <- anis_d1$ycoord_val

    anis_d2 <- transform_anis(
      data = d2,
      data_object$xcoord, data_object$ycoord,
      spcov_params[["rotate"]], spcov_params[["scale"]]
    )

    d2[[data_object$xcoord]] <- anis_d2$xcoord_val
    d2[[data_object$ycoord]] <- anis_d2$ycoord_val
  }

  dist_d1d2_cross <- spdist_vectors(
    data = d1, data2 = d2,
    data_object$xcoord, data_object$ycoord, dim_coords = data_object$dim_coords
  )

  # random effects
  if (!is.null(data_object$randcov_list)) {
    Zs_cross <- lapply(randcov_names, function(x) {
      list(
        ZZt = Matrix::tcrossprod(
          data_object$randcov_list[[d1_index]][[x]]$Z,
          data_object$randcov_list[[d2_index]][[x]]$Z
        )
      )
    })
    names(Zs_cross) <- randcov_names
  } else {
    randcov_params <- NULL
    Zs_cross <- NULL
  }


  # partition matrix
  if (!is.null(data_object$partition_factor)) {
    # finding the formula
    partition_formula <- reformulate(labels(terms(data_object$partition_factor)), intercept = FALSE)
    #
    d1_partition <- Matrix::Matrix(model.matrix(partition_formula, d1), sparse = TRUE)
    d2_partition <- Matrix::Matrix(model.matrix(partition_formula, d2), sparse = TRUE)
    partition_matrix_cross_val <- Matrix::tcrossprod(d1_partition, d2_partition)
  } else {
    partition_matrix_cross_val <- NULL
  }

  cov_d1d2_cross <- cov_matrix_cross(spcov_params, dist_d1d2_cross, randcov_params, Zs_cross, partition_matrix_cross_val)

  part1 <- cholprods_list[[d1_index]]$SqrtSigInv_X
  part2 <- forwardsolve(cholprods_list[[d1_index]]$Sig_lowchol, cov_d1d2_cross)
  part4 <- cholprods_list[[d2_index]]$SqrtSigInv_X
  part3 <- backsolve(t(cholprods_list[[d2_index]]$Sig_lowchol), part4)

  W_ij_half <- Matrix::crossprod(part1, part2) %*% part3
  W_ij <- W_ij_half + t(W_ij_half)
}


get_W_ij_parallel <- function(index_list, cholprods_list, spcov_params,
                              randcov_params, randcov_names, data_object) {
  d1_index <- index_list$d1
  d2_index <- index_list$d2
  get_W_ij(
    d1_index, d2_index, cholprods_list,
    spcov_params, randcov_params,
    randcov_params, data_object
  )
}
