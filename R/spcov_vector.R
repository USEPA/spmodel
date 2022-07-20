#' Create a spatial covariance vector for prediction
#'
#' @param spcov_params An \code{spcov_params} object
#' @param dist_matrix A distance vector specifying the Euclidean distance (splm)
#'   covariances or the neighboring structure (spautor covariances)
#'
#' @return A covariance matrix
#'
#' @noRd
spcov_vector <- function(spcov_params, dist_vector) {
  UseMethod("spcov_vector", spcov_params)
}
########### three parameter geostatistical
# spcov_vector exponential
#' @export
spcov_vector.exponential <- function(spcov_params, dist_vector) {
  spcov_vector_val <- spcov_params[["de"]] * exp(-dist_vector / spcov_params[["range"]])
  spcov_vector_val
}

# spcov_vector spherical
#' @export
spcov_vector.spherical <- function(spcov_params, dist_vector) {
  dist_ratio <- dist_vector / spcov_params[["range"]]
  spcov_vector_val <- spcov_params[["de"]] * (1 - (3 / 2) * dist_ratio + (1 / 2) * dist_ratio^3) * (dist_vector <= spcov_params[["range"]])
  spcov_vector_val
}

# spcov_vector gaussian
#' @export
spcov_vector.gaussian <- function(spcov_params, dist_vector) {
  spcov_vector_val <- spcov_params[["de"]] * exp(-(dist_vector / spcov_params[["range"]])^2)
  spcov_vector_val
}

# spcov_vector triangular
#' @export
spcov_vector.triangular <- function(spcov_params, dist_vector) {
  spcov_vector_val <- spcov_params[["de"]] * (1 - dist_vector / spcov_params[["range"]]) * (dist_vector <= spcov_params[["range"]])
  spcov_vector_val
}

# spcov_vector circular
#' @export
spcov_vector.circular <- function(spcov_params, dist_vector) {
  min_val <- pmin(dist_vector / spcov_params[["range"]], 1)
  spcov_vector_val <- spcov_params[["de"]] * (1 - (2 / pi * (min_val * sqrt(1 - min_val^2) + asin(min_val)))) * (dist_vector <= spcov_params[["range"]])
  spcov_vector_val
}

# spcov_vector none
#' @export
spcov_vector.none <- function(spcov_params, dist_vector) {
  spcov_vector_val <- Matrix(0, nrow = 1, ncol = length(dist_vector), sparse = TRUE) # length dist vector
  spcov_vector_val
}

# spcov_vector cubic
#' @export
spcov_vector.cubic <- function(spcov_params, dist_vector) {
  dist_ratio <- dist_vector / spcov_params[["range"]]
  spcov_vector_val <- spcov_params[["de"]] * (1 - (7 / 1 * dist_ratio^2) + (35 / 4 * dist_ratio^3) - (7 / 2 * dist_ratio^5) + (3 / 4 * dist_ratio^7)) * (dist_vector <= spcov_params[["range"]])
  spcov_vector_val
}

# spcov_vector pentaspherical
#' @export
spcov_vector.pentaspherical <- function(spcov_params, dist_vector) {
  dist_ratio <- dist_vector / spcov_params[["range"]]
  spcov_vector_val <- spcov_params[["de"]] * (1 - (15 / 8 * dist_ratio) + (5 / 4 * dist_ratio^3) - (3 / 8 * dist_ratio^5)) * (dist_vector <= spcov_params[["range"]])
  spcov_vector_val
}

# spcov_vector cosine
#' @export
spcov_vector.cosine <- function(spcov_params, dist_vector) {
  spcov_vector_val <- spcov_params[["de"]] * cos(dist_vector / spcov_params[["range"]])
  spcov_vector_val
}

# spcov_vector wave
#' @export
spcov_vector.wave <- function(spcov_params, dist_vector) {
  spcov_vector_val <- spcov_params[["de"]] * spcov_params[["range"]] * sin(pi * dist_vector / spcov_params[["range"]]) / (pi * dist_vector)
  spcov_vector_val
}

# spcov_vector jbessel
#' @export
spcov_vector.jbessel <- function(spcov_params, dist_vector) {
  dist_product <- dist_vector * spcov_params[["range"]]
  spcov_vector_val <- spcov_params[["de"]] * besselJ(as.matrix(pmin(dist_product, 100000)), 0)
  spcov_vector_val
}

# spcov_vector gravity
#' @export
spcov_vector.gravity <- function(spcov_params, dist_vector) {
  dist_ratio <- dist_vector / spcov_params[["range"]]
  spcov_vector_val <- spcov_params[["de"]] * (1 + dist_ratio^2)^(-1 / 2)
  spcov_vector_val
}

# spcov_vector rquad
#' @export
spcov_vector.rquad <- function(spcov_params, dist_vector) {
  dist_ratio <- dist_vector / spcov_params[["range"]]
  spcov_vector_val <- spcov_params[["de"]] * (1 + dist_ratio^2)^(-1)
  spcov_vector_val
}

# spcov_vector magnetic
#' @export
spcov_vector.magnetic <- function(spcov_params, dist_vector) {
  dist_ratio <- dist_vector / spcov_params[["range"]]
  spcov_vector_val <- spcov_params[["de"]] * (1 + dist_ratio^2)^(-1 / 2)
  spcov_vector_val
}
########### four parameter geostatistical
# spcov_vector matern
#' @export
spcov_vector.matern <- function(spcov_params, dist_vector) {
  eta <- sqrt(2 * spcov_params[["extra"]]) * (dist_vector / spcov_params[["range"]])
  spcov_vector_val <- spcov_params[["de"]] * 2^(1 - spcov_params[["extra"]]) / gamma(spcov_params[["extra"]]) * eta^spcov_params[["extra"]] * besselK(as.matrix(eta), nu = spcov_params[["extra"]])
  spcov_vector_val
}

# spcov_vector cauchy
#' @export
spcov_vector.cauchy <- function(spcov_params, dist_vector) {
  dist_ratio <- dist_vector / spcov_params[["range"]]
  spcov_vector_val <- spcov_params[["de"]] * (1 + dist_ratio^2)^(-spcov_params[["extra"]])
  spcov_vector_val
}

# spcov_vector pexponential
#' @export
spcov_vector.pexponential <- function(spcov_params, dist_vector) {
  spcov_vector_val <- spcov_params[["de"]] * exp(-dist_vector^spcov_params[["extra"]] / spcov_params[["range"]])
  spcov_vector_val
}
