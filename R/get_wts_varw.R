get_wts_varw <- function(family, Xmat, y, w, size, dispersion, cov_lowchol, x0, c0) {

  SigInv <- chol2inv(t(cov_lowchol)) # works on upchol
  SigInv_X <- SigInv %*% Xmat
  cov_betahat <- chol2inv(chol(crossprod(Xmat, SigInv_X)))
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  Ptheta <- SigInv - SigInv_X %*% wts_beta

  d <-  get_d(family, w, y, size, dispersion)
  # and then the gradient vector
  # g <-  d - Ptheta %*% w
  # Next, compute H
  D <- get_D(family, w, y, size, dispersion)
  H <- D - Ptheta
  mHInv <- solve(-H) #chol2inv(chol(Matrix::forceSymmetric(-H))) # solve(-H)

  if (NROW(x0) == 1) {
    wts_pred <- x0 %*% wts_beta + c0 %*% SigInv - (c0 %*% SigInv_X) %*% wts_beta
    var_adj <- as.numeric(wts_pred %*% tcrossprod(mHInv, wts_pred))
  } else {
    var_adj <- vapply(seq_len(NROW(x0)), function(x) {
      x0_new <- x0[x, , drop = FALSE]
      c0_new <- c0[x, , drop = FALSE]
      wts_pred <- x0_new %*% wts_beta + c0_new %*% SigInv - (c0_new %*% SigInv_X) %*% wts_beta
      as.numeric(wts_pred %*% tcrossprod(mHInv, wts_pred))
    }, numeric(1))
  }
  var_adj
}




# get_wts_varw_old <- function(object, newdata, newdata_model, local) {
#
#
#   if (is.null(object$local_index)) {
#     local_index <- rep(1, object$n)
#   } else {
#     local_index <- local_index
#   }
#
#
#   X_list <- split.data.frame(model.matrix(object), local_index)
#   X <- do.call("rbind", X_list)
#   y_list <- split.data.frame(as.matrix(object$y, ncol = 1), local_index)
#   y <- do.call("rbind", y_list)
#   w_list <- split.data.frame(as.matrix(object$w, ncol = 1), local_index)
#   w <- do.call("rbind", w_list)
#   if (is.null(object$size)) {
#     size_list <- NULL
#     size <- NULL
#   } else {
#     size_list <- split.data.frame(as.matrix(object$size, ncol = 1), local_index)
#     size <- do.call("rbind", size_list)
#   }
#
#
#   cov_matrix <- covmatrix(object)
#   cov_matrix_list <- lapply(names(X_list), function(x) {
#     local_row <- which(local_index == x) # numeric to cont conversion happens
#     cov_matrix[local_row, local_row, drop = FALSE]
#   })
#
#   # cholesky products
#   if (local$parallel) {
#     cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
#       cluster_list_element <- list(
#         c = cov_matrix_list[[l]],
#         x = X_list[[l]],
#         y = y_list[[l]]
#       )
#     })
#     cholprods_list <- parallel::parLapply(object$cl, cluster_list, get_cholprods_glm_parallel)
#     names(cholprods_list) <- names(cov_matrix_list)
#   } else {
#     cholprods_list <- mapply(
#       c = cov_matrix_list, x = X_list, y = y_list,
#       function(c, x, y) get_cholprods_glm(c, x, y),
#       SIMPLIFY = FALSE
#     )
#   }
#
#   family <- object$family
#   SigInv_list <- lapply(cholprods_list, function(x) x$SigInv)
#   SigInv <- Matrix::bdiag(SigInv_list)
#   SigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SigInv_X))
#   cov_betahat_Inv <- crossprod(X, SigInv_X)
#   cov_betahat <- chol2inv(chol(cov_betahat_Inv)) # need to compute thi
#   # because vcov() returns adjusted
#   Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
#
#   # find dispersion
#   dispersion <- as.vector(coef(object, type = "dispersion")) # take class away
#
#   d <-  get_d(family, w, y, size, dispersion)
#   # and then the gradient vector
#   g <-  d - Ptheta %*% w
#   # Next, compute H
#   D <- get_D(family, w, y, size, dispersion)
#   D_diag <- diag(D)
#   D_list <- lapply(split(D_diag, local_index), function(x) Diagonal(x = x))
#   # cholesky products
#   if (local$parallel) {
#     cluster_list <- lapply(seq_along(D_list), function(l) {
#       cluster_list_element <- list(
#         D = D_list[[l]],
#         S = SigInv_list[[l]]
#       )
#     })
#     DSigInv_list <- parallel::parLapply(local$cl, cluster_list, get_DSigInv_parallel)
#     names(DSigInv_list) <- names(D_list)
#     DSigInv_Inv_list <- parallel::parLapply(local$cl, DSigInv_list, solve)
#   } else {
#     DSigInv_list <- mapply(
#       D = D_list, S = SigInv_list,
#       function(D, S) get_DSigInv(D, S),
#       SIMPLIFY = FALSE
#     )
#     DSigInv_Inv_list <- lapply(DSigInv_list, solve)
#   }
#
#   DSigInv_Inv <- Matrix::bdiag(DSigInv_Inv_list)
#   mHInv <- -smw_HInv(AInv = DSigInv_Inv, U = SigInv_X, CInv = cov_betahat_Inv)
#   order_val <- unlist(split(seq_len(object$n), local_index), use.names = FALSE)
#   wts_beta <- tcrossprod(cov_betahat, SigInv_X)
#
#   mHInv <- mHInv[order(order_val), order(order_val), drop = FALSE]
#   SigInv <- SigInv[order(order_val), order(order_val), drop = FALSE]
#   SigInv_X <- SigInv_X[order(order_val), , drop = FALSE]
#   wts_beta <- wts_beta[order(order_val)]
#
#
#   c0 <- covmatrix(object, newdata)
#   var_adj <- vapply(seq_len(NROW(newdata)), function(x) {
#     newdata_model_new <- newdata_model[x, , drop = FALSE]
#     c0_new <- c0[x, , drop = FALSE]
#     wts_pred <- newdata_model_new %*% wts_beta + c0_new %*% SigInv - (c0_new %*% SigInv_X) %*% wts_beta
#     var_adj <- as.numeric(wts_pred %*% tcrossprod(mHInv, wts_pred))
#   }, numeric(1))
#
#
#   var_adj
# }
