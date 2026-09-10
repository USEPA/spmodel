#' @noRd
dSig_dtheta_spcov <- function(spcov_params_val, dist_matrix, ...) {
  UseMethod("dSig_dtheta_spcov", spcov_params_val)
}

#' @noRd
dSig_dtheta_spcov.default <- function(spcov_params_val, dist_matrix, ...) {
  stop(
    "Closed-form covariance parameter derivatives are not implemented for spcov_type = \"",
    class(spcov_params_val)[1], "\" yet; use method = \"numeric\" instead.",
    call. = FALSE
  )
}

#' @noRd
# exponential covariance: Sigma = de * R + ie * I, where R_ij = exp(-dist_ij /
# range) is the pure correlation matrix. Each entry below is the closed-form partial
# derivative of that expression with respect to one parameter, holding the
# others fixed:
#  - d(de*R + ie*I)/d(de) = R
#  - d(de*R + ie*I)/d(ie) = I
#  - d(de*R)/d(range) = de * R * (dist/range^2), from the chain rule applied
#    to d/d(range)[exp(-dist/range)] = exp(-dist/range) * (dist/range^2)
dSig_dtheta_spcov.exponential <- function(spcov_params_val, dist_matrix, ...) {
  de <- spcov_params_val[["de"]]
  range <- spcov_params_val[["range"]]
  R <- exp(-dist_matrix / range)
  list(
    de = R,
    ie = Diagonal(nrow(dist_matrix)),
    range = de * R * (dist_matrix / range^2)
  )
}

#' @noRd
# gaussian covariance: Sigma = de * R + ie * I, where R_ij = exp(-(dist_ij /
# range)^2). Same de/ie logic as exponential above; the range derivative
# comes from the chain rule applied to the squared ratio:
#  - d(de*R)/d(range) = de * R * (2*dist^2/range^3), since
#    d/d(range)[exp(-dist^2/range^2)] = exp(-dist^2/range^2) * (2*dist^2/range^3)
dSig_dtheta_spcov.gaussian <- function(spcov_params_val, dist_matrix, ...) {
  de <- spcov_params_val[["de"]]
  range <- spcov_params_val[["range"]]
  R <- exp(-(dist_matrix / range)^2)
  list(
    de = R,
    ie = Diagonal(nrow(dist_matrix)),
    range = de * R * (2 * dist_matrix^2 / range^3)
  )
}

#' @noRd
# spherical covariance: Sigma = de * R + ie * I, where R_ij = (1 - 1.5*ratio +
# 0.5*ratio^3) for ratio = dist_ij/range, truncated to 0 once dist_ij exceeds
# range. de/ie derivatives are the same idea as exponential/gaussian above
# (the truncation indicator is treated as a constant with the value
# one (h <= range) or 0 (h > range); the range
# derivative is then:
#  - d(de*R)/d(range) = de * 1.5 * (dist/range^2) * (1 - ratio^2), from
#    d/d(range)[1 - 1.5*dist/range + 0.5*dist^3/range^3]
#    = 1.5*dist/range^2 - 1.5*dist^3/range^4 = 1.5*(dist/range^2)*(1 - ratio^2)
dSig_dtheta_spcov.spherical <- function(spcov_params_val, dist_matrix, ...) {
  de <- spcov_params_val[["de"]]
  range <- spcov_params_val[["range"]]
  dist_ratio <- dist_matrix / range
  within_range <- dist_matrix <= range
  R <- (1 - 1.5 * dist_ratio + 0.5 * dist_ratio^3) * within_range
  list(
    de = R,
    ie = Diagonal(nrow(dist_matrix)),
    range = de * 1.5 * (dist_matrix / range^2) * (1 - dist_ratio^2) * within_range
  )
}

#' @noRd
# "none"/"ie" covariance: Sigma = ie * I, with de fixed at 0 and range fixed
# at Inf (see spcov_initial_NA()) -- neither actually appears in Sigma, so
#  - d(Sigma)/d(de) = 0
#  - d(Sigma)/d(ie) = I, same as every other spcov_type's ie derivative
#  - d(Sigma)/d(range) = 0
# n is taken from `...` (get_dSig_dtheta_cov() passes it as
# n = context$data_object$n) rather than nrow(dist_matrix), because
# get_cov_gradients_context.splm() leaves dist_matrix NULL for spcov_type
# "none"/"ie" without random effects -- it is genuinely unneeded for the math
# above, so no real distance matrix is built just to size these matrices
dSig_dtheta_spcov.none <- function(spcov_params_val, dist_matrix, ..., n = nrow(dist_matrix)) {
  zero_mat <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  list(
    de = zero_mat,
    ie = Diagonal(n),
    range = zero_mat
  )
}

#' @noRd
# "ie" is a spcov_params()/spcov_initial() alias for "none" (identical
# Sigma = ie * I structure, just a different class label) when 
# splm() or spautor() models are fit (it does differ for spglm() and spgautor(), but
# those functions are not currently used for satterthwaite)
dSig_dtheta_spcov.ie <- dSig_dtheta_spcov.none

#' @noRd
# a random effect contributes sigma_k^2 * Z_k Z_k' to Sigma, linearly in its
# own variance component sigma_k^2 -- so unlike the spatial covariance
# parameters above, every random effect's derivative is simply its own
# (fixed, not parameter-dependent) design product ZZt, already precomputed
# and stored on the data object at fit time
dSig_dtheta_randcov <- function(randcov_Zs, randcov_names) {
  val <- lapply(randcov_names, function(nm) as.matrix(randcov_Zs[[nm]][["ZZt"]]))
  names(val) <- randcov_names
  val
}

# assembles the full list of dSigma/dtheta_k matrices, one per free
# covariance/random-effect parameter, in the same order as
# context$cov_names_free -- get_vcov_theta()'s "closed" branch and
# get_grad_g()'s "closed" branch both use this list directly
get_dSig_dtheta_cov <- function(context, object) {
  UseMethod("get_dSig_dtheta_cov", object)
}

#' @noRd
#' @exportS3Method
get_dSig_dtheta_cov.splm <- function(context, object) {

  randcov_Zs <- if (is.null(context$data_object$randcov_list)) NULL else context$data_object$randcov_list[[1]]

  dSig_spcov <- dSig_dtheta_spcov(context$spcov_params, context$dist_matrix, n = context$data_object$n)[context$spcov_names_free]
  dSig_randcov <- if (!is.null(context$randcov_names_free)) {
    dSig_dtheta_randcov(randcov_Zs, context$randcov_names_free)
  } else {
    list()
  }
  c(dSig_spcov, dSig_randcov)
}

#' @noRd
#' @exportS3Method
get_dSig_dtheta_cov.spautor <- function(context, object) {

  dSig_spcov <- dSig_dtheta_spcov(context$spcov_params, context$data_object$W, n = context$data_object$n)[context$spcov_names_free]
  dSig_randcov <- if (!is.null(context$randcov_names_free)) {
    dSig_dtheta_randcov(context$data_object$randcov_Zs, context$randcov_names_free)
  } else {
    list()
  }
  c(dSig_spcov, dSig_randcov)
}
