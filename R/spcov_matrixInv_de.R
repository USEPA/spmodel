#' Create an covariance matrix inverse for the dependent error portion of
#' CAR and SAR models
#'
#' @param spcov_params An \code{spcov_params} object
#' @param dist_matrix A weights (distance) matrix specifying the data's
#'   neighbor structure
#'
#' @return A covariance matrix
#'
#' @noRd
spcov_matrixInv_de <- function(spcov_params, dist_matrix, ...) {
  UseMethod("spcov_matrixInv_de", spcov_params)
}

# spcov_matrixInv_de car
#' @export
spcov_matrixInv_de.car <- function(spcov_params, dist_matrix, M, ...) {
  if (spcov_params[["de"]] < 0.01) {
    spcov_params[["de"]] <- 0.01
  } # done for invertibility issues arising with cov_initial_search

  # compute row sums of W
  dist_matrix_rowsums <- rowSums(dist_matrix)

  index <- dist_matrix_rowsums == 0
  C_index <- which(!index)
  U_index <- which(index)

  if (length(U_index) > 0) {
    dist_matrix_CC <- dist_matrix[C_index, C_index, drop = FALSE]
    # the off diagonals are not needed because they are assumed zero later
    dist_matrix_UU <- dist_matrix[U_index, U_index, drop = FALSE]

    # compute - rho * W
    cor_matrixInv_CC <- -spcov_params[["range"]] * dist_matrix_CC
    # compute (I - rho * W)
    diag(cor_matrixInv_CC) <- 1 + diag(cor_matrixInv_CC)
    # compute M^{-1} %*% (I - rho * W) (same as R (I - rho * W) * 1 / M)
    cor_matrixInv_CC <- cor_matrixInv_CC * 1 / M[C_index]
    spcov_matrixInv_CC <- 1 / spcov_params[["de"]] * cor_matrixInv_CC

    cor_matrixInv_UU <- Matrix(diag(length(U_index)), sparse = TRUE)
    spcov_matrixInv_UU <- 1 / spcov_params[["extra"]] * cor_matrixInv_UU # no M multiplication here

    spcov_matrixInv_de <- Matrix::Matrix(0, nrow = length(index), ncol = length(index), sparse = TRUE)
    spcov_matrixInv_de[C_index, C_index] <- spcov_matrixInv_CC
    spcov_matrixInv_de[U_index, U_index] <- spcov_matrixInv_UU
  } else {

    # compute - rho * W
    cor_matrixInv_de <- -spcov_params[["range"]] * dist_matrix
    # compute (I - rho * W)
    diag(cor_matrixInv_de) <- 1 + diag(cor_matrixInv_de)
    # compute M^{-1} %*% (I - rho * W) (same as R (I - rho * W) * 1 / M)
    cor_matrixInv_de <- cor_matrixInv_de * 1 / M
    # else do nothing as M is identity
    spcov_matrixInv_de <- 1 / spcov_params[["de"]] * cor_matrixInv_de
  }
  spcov_matrixInv_de
}

# spcov_matrixInv_de sar
#' @export
spcov_matrixInv_de.sar <- function(spcov_params, dist_matrix, M, ...) { # M is not used
  if (spcov_params[["de"]] < 0.01) {
    spcov_params[["de"]] < 0.01
  } # done for invertibility issues arising with cov_initial_search

  # compute row sums of W
  dist_matrix_rowsums <- rowSums(dist_matrix)

  # find connected vs unconnected indices
  index <- dist_matrix_rowsums == 0
  C_index <- which(!index)
  U_index <- which(index)

  if (length(U_index) > 0) {
    dist_matrix_CC <- dist_matrix[C_index, C_index, drop = FALSE]
    # the off diagonals are not needed because they are assumed zero later
    dist_matrix_UU <- dist_matrix[U_index, U_index, drop = FALSE]

    # compute - rho * W
    cor_matrixInv_CC_left <- -spcov_params[["range"]] * dist_matrix_CC
    # compute (I - rho * W)^{-1}
    diag(cor_matrixInv_CC_left) <- 1 + diag(cor_matrixInv_CC_left)
    # compute (I - rho * W) (I - rho * W)'
    cor_matrixInv_CC <- tcrossprod(cor_matrixInv_CC_left, cor_matrixInv_CC_left)
    spcov_matrixInv_CC <- 1 / spcov_params[["de"]] * cor_matrixInv_CC

    cor_matrixInv_UU <- Matrix(diag(length(U_index)), sparse = TRUE)
    spcov_matrixInv_UU <- 1 / spcov_params[["extra"]] * cor_matrixInv_UU

    spcov_matrixInv_de <- Matrix::Matrix(0, nrow = length(index), ncol = length(index), sparse = TRUE)
    spcov_matrixInv_de[C_index, C_index] <- spcov_matrixInv_CC
    spcov_matrixInv_de[U_index, U_index] <- spcov_matrixInv_UU
  } else {

    # compute - rho * W
    cor_matrixInv_de_left <- -spcov_params[["range"]] * dist_matrix
    # compute (I - rho * W)^{-1}
    diag(cor_matrixInv_de_left) <- 1 + diag(cor_matrixInv_de_left)
    # compute (I - rho * W) (I - rho * W)'
    cor_matrixInv_de <- tcrossprod(cor_matrixInv_de_left, cor_matrixInv_de_left)
    spcov_matrixInv_de <- 1 / spcov_params[["de"]] * cor_matrixInv_de
  }

  spcov_matrixInv_de
}
