#' Simulate newdata sequentially, conditional on observed data (Vecchia
#' method), adjusted for latent-process (\code{w}) estimation uncertainty
#' (\code{spglm} models)
#'
#' GLM analog of \code{\link{get_conditional_vecchia}()}, used by
#' \code{conditional.spglm()}. As in \code{get_conditional_vecchia()}, every
#' \code{newdata} location is drawn one at a time, conditional on all observed
#' data plus every earlier-drawn \code{newdata} location. The link-scale
#' latent process \code{w} is held fixed at its fitted value and its own Laplace-approximate posterior uncertainty is
#' instead supplied analytically via \code{var_adj}, exactly as in
#' \code{\link{get_conditional_new_from_base_adjust_glm}()}.
#'
#' \code{var_adj} is not spatially local as every \code{newdata} location's
#' prediction weights span all of the observed data, because the Laplace
#' posterior for \code{w} is a single joint (not spatially-truncatable)
#' distribution over every observed location (see \code{conditional.spglm()}'s
#' one-time factorization of \code{cov_lowchol_mH}). Its contribution to every
#' \code{newdata} location's own (diagonal) predictive variance is therefore
#' computed once for all of \code{newdata} up front (\code{var_adj_diag}),
#' rather than inside the sequential loop. Further investigation should be
#' made into a spatially local version.
#'
#' \code{var_adj} also induces covariance \emph{between} \code{newdata}
#' predictions, since they share the same uncertain \code{w}. Exactly as
#' \code{\link{get_conditional_new_from_base_adjust_glm}()} adds \code{var_adj}
#' only to the newdata-newdata conditional covariance (never to the
#' observed-observed or observed-newdata covariance, since the observed
#' data/base sample is held fixed, not predicted), this cross-newdata
#' covariance is added here only between pairs of pool members that are
#' \emph{both} earlier-drawn \code{newdata} locations (tracked via
#' \code{pool_is_obs}, exactly as in \code{get_conditional_vecchia()}) and never
#' between an observed pool member and anything else. Because this term is
#' folded into the same \code{cov_target_pool}/\code{cov_pool_pool} objects
#' used for the spatial conditioning, it is truncated by the same
#' \code{local_list$method}/\code{size} neighbor selection as the spatial
#' part and there is no separate approximation layer. This also means that,
#' exactly as in \code{get_conditional_vecchia()}, this is mathematically
#' exact when \code{local_list$method == "all"}: sequential conditioning on
#' the "effective" covariance (ordinary spatial covariance, plus \code{var_adj}
#' between newdata pairs) reproduces the same joint distribution as
#' \code{conditional.spglm()}'s \code{"low-rank"} path with
#' \code{method_new = "all"} (i.e. \code{local = FALSE}), which adds the same
#' \code{var_adj} term to its own (single, unblocked) newdata-newdata
#' conditional covariance before a single joint draw.
#'
#' @param object A fitted \code{spglm} model object.
#' @param newdata The (already processed, plain data frame) \code{newdata}
#'   locations to simulate, matching \code{object$xcoord}/\code{object$ycoord}.
#' @param newdata_model The \code{newdata} design matrix, in \code{newdata}'s
#'   original row order.
#' @param base_val A matrix of the (fixed, not simulated) observed-data
#'   link-scale residuals \code{w - X \%*\% new_betahat}, one row per row of
#'   \code{object$obdata} and one column per simulated beta draw (see
#'   \code{conditional.spglm()}).
#' @param local_list The resolved \code{local} list from
#'   \code{\link{get_local_list_conditional}()} (\code{approximation == "vecchia"}).
#' @param samples The number of simulations (columns of \code{base_val}).
#' @param SigInv The precision matrix of the observed data's covariance
#'   matrix (computed over \strong{all} observed data as vecchia never
#'   subsamples it).
#' @param SigInv_X \code{SigInv \%*\% X}, for the full observed design matrix.
#' @param wts_beta \code{cov_betahat \%*\% t(SigInv_X)}, prediction weights for
#'   the fixed effect contribution.
#' @param cov_lowchol_mH The lower triangular Cholesky factor of the negative
#'   Hessian of the joint log-likelihood for \code{w}, factored once over all
#'   observed data (see \code{conditional.spglm()}).
#'
#' @return A matrix of simulated link-scale residuals at \code{newdata}, in
#'   \code{newdata}'s original row order, one column per simulation.
#'
#' @noRd
get_conditional_vecchia_glm <- function(object, newdata, newdata_model, base_val, local_list, samples,
                                         SigInv, SigInv_X, wts_beta, cov_lowchol_mH) {

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

  # var_adj prediction weights for every newdata location at once (not
  # spatially truncated), matching the formula used in
  # get_conditional_new_from_base_adjust_glm() for a single (unsplit) block
  cov_base_new <- covmatrix(object, newdata, cov_type = "obs.pred") # n_obs x n_new
  c0_all <- t(cov_base_new) # n_new x n_obs
  wts_pred_all <- newdata_model %*% wts_beta + c0_all %*% SigInv - (c0_all %*% SigInv_X) %*% wts_beta
  wts_pred_all <- t(wts_pred_all) # n_obs x n_new
  SqrtmHInv_wts_pred_all <- forwardsolve(cov_lowchol_mH, wts_pred_all) # n_obs x n_new
  # var_adj_diag[i] is this newdata location's own var_adj contribution
  # (var_adj_full[i, i], without ever forming the full n_new x n_new matrix)
  var_adj_diag <- colSums(SqrtmHInv_wts_pred_all^2)

  xo <- xcoord_new[ord]
  yo <- ycoord_new[ord]

  Y_ordered <- matrix(NA_real_, n_new, samples)
  Z <- matrix(rnorm(n_new * samples), n_new, samples)

  for (k in seq_len(n_new)) {
    n_new_pool <- k - 1
    pool_x <- if (n_new_pool > 0) c(xcoord_obs, xo[seq_len(n_new_pool)]) else xcoord_obs
    pool_y <- if (n_new_pool > 0) c(ycoord_obs, yo[seq_len(n_new_pool)]) else ycoord_obs
    pool_val <- if (n_new_pool > 0) rbind(base_val, Y_ordered[seq_len(n_new_pool), , drop = FALSE]) else base_val
    pool_is_obs <- c(rep(TRUE, n_obs), rep(FALSE, n_new_pool))
    # for pool members that are earlier-drawn newdata rows, pool_idx holds
    # their (unordered) newdata row index and this lines up directly with
    # SqrtmHInv_wts_pred_all's/var_adj_diag's column indexing below
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
      pool_val <- pool_val[keep, , drop = FALSE]
      dist_target_pool <- dist_target_pool[keep]
      pool_is_obs <- pool_is_obs[keep]
      pool_idx <- pool_idx[keep]
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

    # fold var_adj into the spatial covariance above, but only between pairs
    # that are BOTH "new"/predicted locations. Observed pool
    # members are fixed data, not predictions, so they carry no var_adj, and
    # this is truncated by the same neighbor selection as the spatial part.
    new_pos <- which(!pool_is_obs)
    if (length(new_pos) > 0) {
      new_cols <- pool_idx[new_pos]
      Sqrt_target <- SqrtmHInv_wts_pred_all[, ord[k], drop = FALSE]
      Sqrt_new <- SqrtmHInv_wts_pred_all[, new_cols, drop = FALSE]
      cov_target_pool[new_pos] <- cov_target_pool[new_pos] + as.numeric(crossprod(Sqrt_target, Sqrt_new))
      cov_pool_pool[new_pos, new_pos] <- cov_pool_pool[new_pos, new_pos] + crossprod(Sqrt_new, Sqrt_new)
    }
    target_var <- target_var + var_adj_diag[ord[k]]

    chol_pool <- chol(cov_pool_pool)
    w <- backsolve(chol_pool, forwardsolve(t(chol_pool), cov_target_pool))
    cond_var <- max(target_var - sum(w * cov_target_pool), 0)

    cond_mean <- as.numeric(crossprod(w, pool_val))

    Y_ordered[k, ] <- cond_mean + sqrt(cond_var) * Z[k, ]
  }

  Y <- matrix(NA_real_, n_new, samples)
  Y[ord, ] <- Y_ordered
  Y
}
