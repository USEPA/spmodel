get_fitted_null <- function(w, data_object) {
  fitted_link <- as.numeric(w)
  fitted_response <- invlink(fitted_link, data_object$family, data_object$size)
}

get_fitted_spglm <- function(w_list, betahat, spcov_params, data_object, eigenprods_list,
                           dist_matrix_list, randcov_params = NULL) {

  fitted_link <- unname(do.call("c", w_list)) # unlist(w_list, use.names = FALSE)
  fitted_response <- invlink(fitted_link, data_object$family, data_object$size)

  SigInv_r_list <- mapply(x = eigenprods_list, w = w_list, function(x, w) x$SigInv %*% as.matrix(w, ncol = 1) - x$SigInv_X %*% betahat,
                          SIMPLIFY = FALSE)

  # cov params no de   (set ie portion to zero because BLUP only uses cov(dependent error))
  spcov_params_de_only <- spcov_params
  spcov_params_de_only[["ie"]] <- 0
  spcov_matrix_de_only_list <- lapply(
    dist_matrix_list,
    function(x) spcov_matrix(spcov_params = spcov_params_de_only, dist_matrix = x)
  )


  ### cov(dependent error) * Z (identity) * siginv (y - x beta)
  fitted_de <- as.numeric(do.call("rbind", mapply(
    s = spcov_matrix_de_only_list, r = SigInv_r_list,
    function(s, r) s %*% r, SIMPLIFY = FALSE
  )))

  ### cov(independent error) is zero so this gives
  ### sigma^2(independent) * Identity * Z (identity) * siginv (y - x beta)
  fitted_ie <- as.numeric(spcov_params[["ie"]] * do.call("rbind", SigInv_r_list))

  ## fitted random effects
  if (is.null(names(randcov_params))) {
    fitted_randcov <- NULL
  } else {
    fitted_randcov <- lapply(names(randcov_params), function(x) {
      fitted_val <- randcov_params[[x]] * do.call("rbind", mapply(
        z = data_object$randcov_list,
        r = SigInv_r_list,
        function(z, r) {
          crossprod(z[[x]][["Z"]], r)
        }
      ))
      fitted_val <- tapply(fitted_val, rownames(fitted_val), function(x) {
        val <- mean(x[x != 0])
        if (length(val) == 0) { # replace if all zeros somehow
          val <- rep(0, length(x))
          names(val) <- names(x)
        }
        val
      })
      # all combinations yields values with many zeros -- don't want to include these in the mean
      names_fitted_val <- rownames(fitted_val)
      fitted_val <- as.numeric(fitted_val)
      names(fitted_val) <- names_fitted_val
      fitted_val
    })
    names(fitted_randcov) <- names(randcov_params)
  }

  # call latent link?
  fitted_values <- list(
    response = fitted_response,
    link = fitted_link,
    spcov = list(de = fitted_de, ie = fitted_ie),
    randcov = fitted_randcov
  )

}

get_fitted_spgautor <- function(w, betahat, spcov_params, data_object, eigenprods,
                           dist_matrix_list, randcov_params = NULL) {

  fitted_link <- as.numeric(w)
  fitted_response <- invlink(fitted_link, data_object$family, data_object$size)

  dist_matrix <- data_object$W[data_object$observed_index, data_object$observed_index, drop = FALSE]
  M <- data_object$M[data_object$observed_index]

  SigInv_r <- eigenprods$SigInv %*% w - eigenprods$SigInv_X %*% betahat

  # cov params no de   (set ie portion to zero because BLUP only uses cov(dependent error))
  spcov_params_de_only <- spcov_params
  spcov_params_de_only[["ie"]] <- 0
  spcov_matrix_de_only <- spcov_matrix(spcov_params = spcov_params_de_only, dist_matrix = dist_matrix, M = M)

  if (!is.null(data_object$partition_factor)) {
    spcov_matrix_de_only <- spcov_matrix_de_only * data_object$partition_matrix[data_object$observed_index, data_object$observed_index, drop = FALSE]
  }

  ### cov(dependent error) * Z (identity) * siginv (y - x beta)
  fitted_de <- spcov_matrix_de_only %*% SigInv_r

  ### cov(independent error) is zero so this gives
  ### sigma^2(independent) * Identity * Z (identity) * siginv (y - x beta)
  fitted_ie <- spcov_params[["ie"]] * SigInv_r

  ## fitted random effects
  if (is.null(names(randcov_params))) {
    fitted_randcov <- NULL
  } else {
    if (is.null(data_object$partition_factor)) {
      ob_randcov_Zs <- get_randcov_Zs(data_object$obdata, names(randcov_params), ZZt = FALSE)
      fitted_randcov <- lapply(names(randcov_params), function(x) {
        fitted_val <- randcov_params[[x]] * crossprod(ob_randcov_Zs[[x]][["Z"]], SigInv_r)
        names_fitted_val <- rownames(fitted_val)
        fitted_val <- as.vector(fitted_val)
        names(fitted_val) <- names_fitted_val
        fitted_val
      })
      names(fitted_randcov) <- names(randcov_params)
    } else {
      index <- unname(model.response(model.frame(reformulate("1", response = labels(terms(data_object$partition_factor))),
                                                 data = data_object$obdata
      )))
      index_val <- unique(index)
      ob_randcov_Zs <- get_randcov_Zs(data_object$obdata, names(randcov_params), ZZt = FALSE)
      fitted_randcov <- lapply(names(randcov_params), function(x) {
        fitted_val <- lapply(index_val, function(y) {
          row_val <- y == index
          fitted_vals <- randcov_params[[x]] *
            crossprod(ob_randcov_Zs[[x]][["Z"]][row_val, , drop = FALSE], SigInv_r[row_val, , drop = FALSE])
        })
        fitted_val <- do.call("rbind", fitted_val)
        fitted_val <- tapply(fitted_val, rownames(fitted_val), function(x) {
          val <- mean(x[x != 0])
          if (length(val) == 0) { # replace if all zeros somehow
            val <- rep(0, length(x))
            names(val) <- names(x)
          }
          val
        })
        names_fitted_val <- rownames(fitted_val)
        fitted_val <- as.vector(fitted_val)
        names(fitted_val) <- names_fitted_val
        fitted_val
      })
      names(fitted_randcov) <- names(randcov_params)
    }
  }

  # call latent link?
  fitted_values <- list(
    response = fitted_response,
    link = fitted_link,
    spcov = list(de = fitted_de, ie = fitted_ie),
    randcov = fitted_randcov
  )

}

invlink <- function(fitted_link, family, size) {

  if (family == "poisson") {
    fitted <- exp(fitted_link)
  } else if (family == "binomial") {
    fitted <- size * expit(fitted_link)
  } else if (family == "nbinomial") {
    fitted <- exp(fitted_link)
  } else if (family == "Gamma") {
    # fitted <- 1 / fitted_link
    fitted <- exp(fitted_link)
  } else if (family == "inverse.gaussian") {
    fitted <- exp(fitted_link)
  } else if (family == "beta") {
    fitted <- expit(fitted_link)
  }
  fitted
}
