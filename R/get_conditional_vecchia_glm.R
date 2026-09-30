#' Simulate spatial residuals sequentially from coupled latent and beta draws
#'
#' @param object A fitted spglm object.
#' @param newdata Processed prediction data in original row order.
#' @param base_val Sampled observed residuals w - X beta, one draw per column.
#' @param local_list Resolved Vecchia settings.
#' @param samples Number of draws.
#' @return Simulated prediction residuals in original row order.
#' @noRd
get_conditional_vecchia_glm <- function(object, newdata, base_val, local_list, samples) {
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
  total_var <- spcov_target_var(spcov_val, object$diagtol) + sum(coef(object, type = "randcov"))
  has_randstruct <- !is.null(object$random) || !is.null(object$partition_factor) || object$anisotropy
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
    n_new_pool <- k - 1
    pool_n <- n_obs + n_new_pool
    pool_x <- if (n_new_pool > 0) c(xcoord_obs, xo[seq_len(n_new_pool)]) else xcoord_obs
    pool_y <- if (n_new_pool > 0) c(ycoord_obs, yo[seq_len(n_new_pool)]) else ycoord_obs
    pool_is_obs <- c(rep(TRUE, n_obs), rep(FALSE, n_new_pool))
    pool_idx <- c(seq_len(n_obs), if (n_new_pool > 0) ord[seq_len(n_new_pool)])

    dist_target_pool <- as.numeric(spdist_vectors2(xo[k], yo[k], pool_x, pool_y, sparse = FALSE))
    npool <- length(pool_x)

    if (neighbor_method != "all" && npool > size) {
      if (neighbor_method == "distance") {
        keep <- order(dist_target_pool)[seq_len(size)]
      } else {
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

    chol_pool <- chol(cov_pool_pool)
    w <- backsolve(chol_pool, forwardsolve(t(chol_pool), cov_target_pool))
    cond_var <- target_var - sum(w * cov_target_pool)
    if (cond_var < -sqrt(.Machine$double.eps) * abs(target_var)) {
      stop("The conditional spatial generalized linear model spatial variance is negative; check the fitted covariance.", call. = FALSE)
    }
    cond_var <- max(cond_var, 0)

    cond_mean <- as.numeric(crossprod(w, pool_val))

    Y_ordered[k, ] <- cond_mean + sqrt(cond_var) * Z[k, ]
    pool_val_full[n_obs + k, ] <- Y_ordered[k, ]
  }

  Y <- matrix(NA_real_, n_new, samples)
  Y[ord, ] <- Y_ordered
  Y
}
