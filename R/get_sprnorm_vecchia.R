#' Simulate every location sequentially (Vecchia method), unconditionally
#'
#' Draws a spatial Gaussian field at every row of \code{object$obdata} one
#' at a time, in \code{local_list$order}, each conditional on every
#' earlier-drawn location. This
#' is the unconditional analog of \code{\link{get_conditional_vecchia}()}
#' (which additionally conditions on already-known observed data); see that
#' function's documentation for the shared exactness argument (a Cholesky-
#' decomposition-as-sequential-conditioning identity) and the shared
#' per-location-not-per-sample vectorization strategy.
#'
#' @param object A fitted \code{splm} model object (built from known/given
#'   covariance parameters purely as a convenient wrapper around
#'   \code{object$obdata}/\code{object$xcoord}/\code{object$ycoord}/
#'   \code{object$random}/\code{object$partition_factor} and \code{covmatrix()} --
#'   see \code{\link{sprnorm.exponential}()}).
#' @param local_list The resolved \code{local} list from
#'   \code{\link{get_local_list_simulation}()} (\code{approximation == "vecchia"}),
#'   supplying \code{order}, \code{method}, and \code{size}.
#' @param samples The number of simulations.
#'
#' @return A matrix of simulated values, in \code{object$obdata}'s original
#'   row order, one column per simulation.
#'
#' @noRd
get_sprnorm_vecchia <- function(object, local_list, samples) {

  xcoord <- object$xcoord
  ycoord <- object$ycoord
  obdata <- object$obdata

  xcoord_val <- obdata[[xcoord]]
  ycoord_val <- obdata[[ycoord]]

  n <- NROW(obdata)

  ord <- local_list$order
  neighbor_method <- local_list$method
  size <- local_list$size

  spcov_val <- coef(object, type = "spcov")
  total_var <- spcov_val[["de"]] + spcov_val[["ie"]] + sum(coef(object, type = "randcov"))
  has_randstruct <- !is.null(object$random) || !is.null(object$partition_factor)
  if (has_randstruct) {
    keep_cols <- unique(c(xcoord, ycoord, all.vars(object$random), all.vars(object$partition_factor)))
  }

  xo <- xcoord_val[ord]
  yo <- ycoord_val[ord]

  Y_ordered <- matrix(NA_real_, n, samples)
  Z <- matrix(rnorm(n * samples), n, samples)

  # first point in the order has no pool to condition on and is drawn directly
  # from its own marginal variance (covmatrix()-derived when random effects/
  # partition factors are present, since a random slope can make marginal
  # variance location-specific; total_var otherwise)
  if (has_randstruct) {
    target_row_1 <- obdata[ord[1], keep_cols, drop = FALSE]
    own_var_1 <- as.numeric(covmatrix(object, newdata = target_row_1, cov_type = "pred.pred"))
  } else {
    own_var_1 <- total_var
  }
  Y_ordered[1, ] <- sqrt(own_var_1) * Z[1, ]

  for (k in 2:n) {
    # conditioning pool: every earlier-drawn location (in the simulation
    # order) and never subsampled, per the design decision that distinguishes
    # "vecchia" from "low-rank"
    n_pool <- k - 1
    pool_x <- xo[seq_len(n_pool)]
    pool_y <- yo[seq_len(n_pool)]
    pool_val <- Y_ordered[seq_len(n_pool), , drop = FALSE]
    pool_idx <- ord[seq_len(n_pool)]

    dist_target_pool <- as.numeric(spdist_vectors2(xo[k], yo[k], pool_x, pool_y, sparse = FALSE))

    if (neighbor_method != "all" && n_pool > size) {
      if (neighbor_method == "distance") {
        keep <- order(dist_target_pool)[seq_len(size)]
      } else {
        # "covariance": rank by |covariance| rather than raw covariance, for
        # spcov_types with negative covariance lobes (wave, cosine, jbessel)
        # (see the identical reasoning in predict.R/decorrelate_data.R).
        # Spatial-only, even when random effects/partition factors are
        # present (see get_conditional_vecchia()'s equivalent note).
        cov_target_pool_full <- as.numeric(cov_vector(spcov_val, dist_target_pool))
        keep <- order(abs(cov_target_pool_full), decreasing = TRUE)[seq_len(size)]
      }
      pool_x <- pool_x[keep]
      pool_y <- pool_y[keep]
      pool_val <- pool_val[keep, , drop = FALSE]
      dist_target_pool <- dist_target_pool[keep]
      pool_idx <- pool_idx[keep]
    }

    if (has_randstruct) {
      # build the (target, neighbor_1, ..., neighbor_m) data frame in the
      # EXACT order of pool_val so cov_full's rows/columns 2:(m+1) line up
      # positionally with pool_val's rows after truncation
      neighbor_rows <- obdata[pool_idx, keep_cols, drop = FALSE]
      target_row <- obdata[ord[k], keep_cols, drop = FALSE]
      combined <- rbind(target_row, neighbor_rows)
      cov_full <- as.matrix(covmatrix(object, newdata = combined, cov_type = "pred.pred"))
      cov_target_pool <- cov_full[1, -1]
      cov_pool_pool <- cov_full[-1, -1, drop = FALSE]
      target_var <- cov_full[1, 1]
    } else {
      cov_target_pool <- as.numeric(cov_vector(spcov_val, dist_target_pool))
      dist_pool_pool <- spdist_vectors2(pool_x, pool_y, pool_x, pool_y, sparse = FALSE)
      cov_pool_pool <- as.matrix(cov_matrix2(spcov_val, dist_matrix = dist_pool_pool))
      target_var <- total_var
    }

    # conditioning weights depend only on locations/covariance parameters
    # computed once per location, reused across every sample below
    chol_pool <- chol(cov_pool_pool)
    w <- backsolve(chol_pool, forwardsolve(t(chol_pool), cov_target_pool))
    cond_var <- max(target_var - sum(w * cov_target_pool), 0)

    # vectorized across every sample (column) of pool_val at once
    cond_mean <- as.numeric(crossprod(w, pool_val))

    Y_ordered[k, ] <- cond_mean + sqrt(cond_var) * Z[k, ]
  }

  Y <- matrix(NA_real_, n, samples)
  Y[ord, ] <- Y_ordered
  Y
}
