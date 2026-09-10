#' Simulate newdata sequentially, conditional on observed data (Vecchia method)
#'
#' Given a Gaussian process already observed (residualized against a
#' simulated beta draw) at every observed location (\code{base_val}), draws
#' values at \code{newdata} locations one at a time, in \code{local_list$order},
#' each conditional on \strong{all} observed data plus every earlier-drawn
#' \code{newdata} location and not a single shared base sample, and not blocks
#' treated as conditionally independent (contrast \code{\link{get_conditional_new_from_base}()},
#' the \code{"low-rank"} type). This is mathematically exact when
#' \code{local_list$method == "all"} (no truncation): sequential
#' conditioning via the chain rule for Gaussians is the same
#' Cholesky-decomposition identity used to justify \code{"low-rank"}'s own
#' block draws, applied to the (already conditional on
#' observed data) joint distribution of \code{newdata} itself, rather than to
#' an unconditional field. Truncating to the \code{local_list$size} nearest/
#' most-correlated candidates in the conditioning pool (\code{method ==
#' "distance"}/\code{"covariance"}) turns this into the approximate
#' Vecchia method, which the observed-data pool is never subsampled for (per
#' the design decision that \code{"vecchia"} exists specifically to avoid the
#' conditional independence bias that \code{"low-rank"} can have).
#'
#' The per-point conditioning weights (\code{w}) and conditional variance
#' depend only on locations and covariance parameters, not on the simulated
#' values themselves, so they are computed once per \code{newdata} location
#' and applied to every sample (column) of the conditioning pool's values at
#' once via a single matrix product. Only the loop over \code{newdata}
#' locations is genuinely sequential (each location's neighbor pool can
#' include earlier-drawn locations), matching every sample simultaneously at
#' each step rather than looping over samples too.
#'
#' If \code{object} has a \code{random} or \code{partition_factor} structure,
#' the pool covariance is instead built by calling \code{\link{covmatrix}()}
#' on a small ad hoc data frame combining the target location with its
#' (possibly truncated) neighbor pool and reusing already-tested code for the
#' random-effect/partition-factor/anisotropy logic rather than re-deriving it
#' here, at the cost of \code{covmatrix()}'s per-call overhead being paid once
#' per \code{newdata} location instead of amortized. Neighbor \emph{ranking}
#' (when \code{local_list$method == "covariance"}) always uses the cheaper
#' spatial-only covariance, matching the existing "local" neighbor-selection
#' precedent elsewhere in the package (\code{predict()}, \code{decorrelate()}) --
#' only the final small conditioning covariance for the truncated set is
#' random-effect/partition-factor compliant.
#'
#' @param object A fitted \code{splm} model object.
#' @param newdata The (already processed, plain data frame) \code{newdata}
#'   locations to simulate, matching \code{object$xcoord}/\code{object$ycoord}.
#' @param base_val A matrix of simulated (mean-zero) observed-data residuals,
#'   one row per row of \code{object$obdata} (in the same order) and one
#'   column per simulation (see \code{\link{get_conditional_new_from_base}()}).
#' @param local_list The resolved \code{local} list from
#'   \code{\link{get_local_list_conditional}()} (\code{approximation == "vecchia"}),
#'   supplying \code{order} (the \code{newdata} simulation order),
#'   \code{method} (the neighbor-selection rule: \code{"all"}/\code{"distance"}/
#'   \code{"covariance"}, same convention as \code{predict()}'s \code{local$method}),
#'   and \code{size}.
#' @param samples The number of simulations (columns of \code{base_val}).
#'
#' @return A matrix of simulated residuals at \code{newdata}, in
#'   \code{newdata}'s original row order, one column per simulation.
#'
#' @noRd
get_conditional_vecchia <- function(object, newdata, base_val, local_list, samples) {

  xcoord <- object$xcoord
  ycoord <- object$ycoord
  obdata <- object$obdata

  xcoord_obs <- obdata[[xcoord]]
  ycoord_obs <- obdata[[ycoord]]
  xcoord_new <- newdata[[xcoord]]
  ycoord_new <- newdata[[ycoord]]

  n_obs <- NROW(obdata)
  n_new <- NROW(newdata)

  ord <- local_list$order
  neighbor_method <- local_list$method
  size <- local_list$size

  spcov_val <- coef(object, type = "spcov")
  total_var <- spcov_val[["de"]] + spcov_val[["ie"]] + sum(coef(object, type = "randcov"))
  has_randstruct <- !is.null(object$random) || !is.null(object$partition_factor)
  if (has_randstruct) {
    keep_cols <- unique(c(xcoord, ycoord, all.vars(object$random), all.vars(object$partition_factor)))
  }

  xo <- xcoord_new[ord]
  yo <- ycoord_new[ord]

  Y_ordered <- matrix(NA_real_, n_new, samples)
  Z <- matrix(rnorm(n_new * samples), n_new, samples)

  # DO NOT USE RBIND IT IS EXTREMELY SLOW,
  # instead pre-allocate the obs-then-new pool once and write each drawn row
  # in place, so a given iteration's pool is always just a slice/subset of
  # this one matrix.
  pool_val_full <- matrix(NA_real_, n_obs + n_new, samples)
  pool_val_full[seq_len(n_obs), ] <- base_val

  for (k in seq_len(n_new)) {
    # conditioning pool: all observed data plus every earlier-drawn newdata
    # location (in the simulation order) and not subsampled, per the design
    # decision that distinguishes "vecchia" from "low-rank"
    n_new_pool <- k - 1
    pool_n <- n_obs + n_new_pool
    pool_x <- if (n_new_pool > 0) c(xcoord_obs, xo[seq_len(n_new_pool)]) else xcoord_obs
    pool_y <- if (n_new_pool > 0) c(ycoord_obs, yo[seq_len(n_new_pool)]) else ycoord_obs
    # tracks provenance (observed row vs. earlier-drawn newdata row) so the
    # has_randstruct branch below can slice the right source data frame,
    # kept in lockstep with pool_x/pool_y/pool_val through truncation
    pool_is_obs <- c(rep(TRUE, n_obs), rep(FALSE, n_new_pool))
    pool_idx <- c(seq_len(n_obs), if (n_new_pool > 0) ord[seq_len(n_new_pool)])

    dist_target_pool <- as.numeric(spdist_vectors2(xo[k], yo[k], pool_x, pool_y, sparse = FALSE))
    npool <- length(pool_x)

    if (neighbor_method != "all" && npool > size) {
      if (neighbor_method == "distance") {
        keep <- order(dist_target_pool)[seq_len(size)]
      } else {
        # "covariance": rank by |covariance| rather than raw covariance, for
        # spcov_types with negative covariance lobes (wave, cosine, jbessel)
        # (see the identical reasoning in predict.R/decorrelate_data.R).
        # Spatial-only, even when random effects/partition factors are
        # present.
        cov_target_pool_full <- as.numeric(cov_vector(spcov_val, dist_target_pool))
        keep <- order(abs(cov_target_pool_full), decreasing = TRUE)[seq_len(size)]
      }
      pool_x <- pool_x[keep]
      pool_y <- pool_y[keep]
      pool_val <- pool_val_full[keep, , drop = FALSE]
      dist_target_pool <- dist_target_pool[keep]
      pool_is_obs <- pool_is_obs[keep]
      pool_idx <- pool_idx[keep]
    } else {
      pool_val <- pool_val_full[seq_len(pool_n), , drop = FALSE]
    }

    if (has_randstruct) {
      # build the (target, neighbor_1, ..., neighbor_m) data frame in the
      # EXACT order of pool_val (one row per lapply() element, rather than
      # rbinding an obs-block and a new-block separately) so cov_full's rows/
      # columns 2:(m+1) line up positionally with pool_val's rows after
      # truncation scrambles the original obs-then-new order
      neighbor_rows <- do.call(rbind, lapply(seq_along(pool_idx), function(i) {
        if (pool_is_obs[i]) obdata[pool_idx[i], keep_cols, drop = FALSE] else newdata[pool_idx[i], keep_cols, drop = FALSE]
      }))
      target_row <- newdata[ord[k], keep_cols, drop = FALSE]
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
    # computed once per newdata location, reused across every sample below
    chol_pool <- chol(cov_pool_pool)
    w <- backsolve(chol_pool, forwardsolve(t(chol_pool), cov_target_pool))
    cond_var <- max(target_var - sum(w * cov_target_pool), 0)

    # vectorized across every sample (column) of pool_val at once
    cond_mean <- as.numeric(crossprod(w, pool_val))

    Y_ordered[k, ] <- cond_mean + sqrt(cond_var) * Z[k, ]
    pool_val_full[n_obs + k, ] <- Y_ordered[k, ]
  }

  Y <- matrix(NA_real_, n_new, samples)
  Y[ord, ] <- Y_ordered
  Y
}
