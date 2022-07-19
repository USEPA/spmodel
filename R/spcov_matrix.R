#' Create a spatial covariance matrix
#'
#' @param spcov_params An \code{spcov_params} object
#' @param dist_matrix A distance matrix specifying the Euclidean distance (splm)
#'   covariances or the neighboring structure (spautor covariances)
#'
#' @return A covariance matrix
#'
#' @noRd
spcov_matrix <- function(spcov_params, dist_matrix, ...) {
  UseMethod("spcov_matrix", spcov_params)
}
########### three parameter geostatistical

# spcov_matrix exponential
#' @export
spcov_matrix.exponential <- function(spcov_params, dist_matrix) {
  spcov_matrix_val <- spcov_params[["de"]] * exp(-dist_matrix / spcov_params[["range"]])
  # numerical stability for positive definiteness
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix spherical
#' @export
spcov_matrix.spherical <- function(spcov_params, dist_matrix) {
  dist_ratio <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * (1 - (3 / 2) * dist_ratio + (1 / 2) * dist_ratio^3) * (dist_matrix <= spcov_params[["range"]])
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix gaussian
#' @export
spcov_matrix.gaussian <- function(spcov_params, dist_matrix) {
  spcov_matrix_val <- spcov_params[["de"]] * exp(-(dist_matrix / spcov_params[["range"]])^2)
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix triangular
#' @export
spcov_matrix.triangular <- function(spcov_params, dist_matrix) {
  spcov_matrix_val <- spcov_params[["de"]] * (1 - dist_matrix / spcov_params[["range"]]) * (dist_matrix <= spcov_params[["range"]])
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix circular
#' @export
spcov_matrix.circular <- function(spcov_params, dist_matrix) {
  min_val <- pmin(dist_matrix / spcov_params[["range"]], 1) # equivalent to below but computationally simpler -- no NaN
  # min_val <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * (1 - (2 / pi * (min_val * sqrt(1 - min_val^2) + asin(min_val)))) * (dist_matrix <= spcov_params[["range"]])
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix none
#' @export
spcov_matrix.none <- function(spcov_params, dist_matrix) {
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]]) # de zero here
  spcov_matrix_val <- diag(rep(spcov_params[["ie"]]), NROW(dist_matrix))
  spcov_matrix_val
}

# spcov_matrix cubic
#' @export
spcov_matrix.cubic <- function(spcov_params, dist_matrix) {
  dist_ratio <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * (1 - (7 / 1 * dist_ratio^2) + (35 / 4 * dist_ratio^3) - (7 / 2 * dist_ratio^5) + (3 / 4 * dist_ratio^7)) * (dist_matrix <= spcov_params[["range"]])
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix pentaspherical
#' @export
spcov_matrix.pentaspherical <- function(spcov_params, dist_matrix) {
  dist_ratio <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * (1 - (15 / 8 * dist_ratio) + (5 / 4 * dist_ratio^3) - (3 / 8 * dist_ratio^5)) * (dist_matrix <= spcov_params[["range"]])
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix cosine
#' @export
spcov_matrix.cosine <- function(spcov_params, dist_matrix) {
  spcov_matrix_val <- spcov_params[["de"]] * cos(dist_matrix / spcov_params[["range"]])
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix wave
#' @export
spcov_matrix.wave <- function(spcov_params, dist_matrix) {
  dist_ratio <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * sin(dist_ratio) / (dist_ratio)
  dist_matrix_zero <- which(dist_matrix == 0)
  spcov_matrix_val[dist_matrix_zero] <- spcov_params[["de"]]
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- spcov_params[["de"]] + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix jbessel
#' @export
spcov_matrix.jbessel <- function(spcov_params, dist_matrix) {
  # besselJ has problem with besselJ(Inf, 0) or anything larger than (100000, 0)
  # pmin() uses a generic based on the type of the first argument
  # vector in first argument goes to vector but matrix in first argument goes to matrix
  dist_product <- dist_matrix * spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * besselJ(as.matrix(pmin(dist_product, 100000)), 0)
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val # pmin operation may take sparse matrix
}

# spcov_matrix gravity
#' @export
spcov_matrix.gravity <- function(spcov_params, dist_matrix) {
  dist_ratio <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * (1 + dist_ratio^2)^(-1 / 2)
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix rquad
#' @export
spcov_matrix.rquad <- function(spcov_params, dist_matrix) {
  dist_ratio <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * (1 + dist_ratio^2)^(-1)
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix magnetic
#' @export
spcov_matrix.magnetic <- function(spcov_params, dist_matrix) {
  dist_ratio <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * (1 + dist_ratio^2)^(-3 / 2)
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}
########### four parameter geostatistical
# spcov_matrix.matern <- function(spcov_params, dist_matrix) {
#   eta <- sqrt(2 * spcov_params[["extra"]]) * (dist_matrix / spcov_params[["range"]])
#   spcov_matrix_val <- spcov_params[["de"]] * 2^(1 - spcov_params[["extra"]]) / gamma(spcov_params[["extra"]]) * eta^spcov_params[["extra"]] * besselK(eta, nu = spcov_params[["extra"]])
#   # with large smoothness, Bessel becomes infinity and is multiplied by something
#   # very small -- so this is fixed to maximum correlation
#   # also -- the bessel function can be very unstable for eta near zero
#   # in my experience, there are some times where very small * besselK(very small, large)
#   # still yields Inf instead of nan, which implies that the limit as eta goes to zero
#   # yields a covariance that should be nearly maximal for the dependent portion (I think
#   # the limit is related to exp(eta) -- just checked, it is)
#   spcov_matrix_val[is.nan(spcov_matrix_val) | is.infinite(spcov_matrix_val)] <- spcov_params[["de"]]
#   if (spcov_params[["ie"]] < 1e-5) {
#     spcov_params[["ie"]] <- 1e-5
#   }
#   diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
#   spcov_matrix_val
# }

# spcov_matrix matern
#' @export
spcov_matrix.matern <- function(spcov_params, dist_matrix) {
  eta <- sqrt(2 * spcov_params[["extra"]]) * (dist_matrix / spcov_params[["range"]])
  spcov_matrix_val <- spcov_params[["de"]] * 2^(1 - spcov_params[["extra"]]) / gamma(spcov_params[["extra"]]) * eta^spcov_params[["extra"]] * besselK(as.matrix(eta), nu = spcov_params[["extra"]]) # eta as sparse matrix causes error
  dist_matrix_zero <- which(dist_matrix == 0)
  spcov_matrix_val[dist_matrix_zero] <- spcov_params[["de"]]
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- spcov_params[["de"]] + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix cauchy
#' @export
spcov_matrix.cauchy <- function(spcov_params, dist_matrix) {
  dist_ratio <- dist_matrix / spcov_params[["range"]]
  spcov_matrix_val <- spcov_params[["de"]] * (1 + dist_ratio^2)^(-spcov_params[["extra"]])
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix pexponential
#' @export
spcov_matrix.pexponential <- function(spcov_params, dist_matrix) {
  spcov_matrix_val <- spcov_params[["de"]] * exp(-dist_matrix^spcov_params[["extra"]] / spcov_params[["range"]])
  spcov_params[["ie"]] <- max(spcov_params[["ie"]], 1e-4 * spcov_params[["de"]])
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}
########### autoregressive

# spcov_matrix car
#' @export
spcov_matrix.car <- function(spcov_params, dist_matrix, M) {
  # find inverse of dependent error (1 / sigma^2 * (I - rho W))
  SigInv_de_val <- spcov_matrixInv_de.car(spcov_params, dist_matrix, M)
  # find covariance of dependent error (sigma^2 * (I - rho W)^-1
  SigInv_de_val_upchol <- chol(forceSymmetric(SigInv_de_val))
  spcov_matrix_val <- chol2inv(SigInv_de_val_upchol)
  ## this approach seems faster (not tested) than SMW to find inverse of
  ## de + ie like is done with spautor_cov_MatrixInv
  # adding independent error variance
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}

# spcov_matrix sar
#' @export
spcov_matrix.sar <- function(spcov_params, dist_matrix, M) { # M is not used
  # find inverse of dependent error (1 / sigma^2 * (I - rho W))
  SigInv_de_val <- spcov_matrixInv_de.sar(spcov_params, dist_matrix)
  # find covariance of dependent error (sigma^2 * (I - rho W)^-1
  SigInv_de_val_upchol <- chol(SigInv_de_val)
  spcov_matrix_val <- chol2inv(SigInv_de_val_upchol)
  # adding independent error variance
  diag(spcov_matrix_val) <- diag(spcov_matrix_val) + spcov_params[["ie"]]
  spcov_matrix_val
}
