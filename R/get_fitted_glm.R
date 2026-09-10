#' Get response-scale fitted values from the latent \code{w} vector alone
#'
#' @param w The latent (link-scale) predictor vector
#' @param data_object The data object
#'
#' @return Response-scale fitted values
#'
#' @noRd
get_fitted_null <- function(w, data_object) {
  # used for null/intercept-only-style fits: w already is the link-scale
  # predictor, so just add the offset and invert the link
  fitted_link <- as.numeric(w)
  if (!is.null(data_object$offset)) {
    fitted_link <- fitted_link + data_object$offset
  }
  invlink(fitted_link, data_object$family, data_object$size)
}

#' Get fitted values for an \code{spglm()} model
#'
#' @param w_list A list of latent (link-scale) predictor vectors, one per partition
#' @param betahat Vector of fixed effects
#' @param spcov_params A \code{spcov_params} object
#' @param data_object The data object
#' @param eigenprods_list A list of cholesky products for each grouping
#' @param dist_matrix_list A list of distance matrices
#' @param randcov_params A \code{randcov_params} object
#'
#' @return A list of fitted values
#'
#' @noRd
get_fitted_spglm <- function(w_list, betahat, spcov_params, data_object, eigenprods_list,
                             dist_matrix_list, randcov_params = NULL) {
  fitted_link <- unname(do.call("c", w_list)) # unlist(w_list, use.names = FALSE)
  # add offset
  if (!is.null(data_object$offset)) {
    fitted_link <- fitted_link + as.vector(data_object$offset)
  }
  # invert the link function to get fitted values on the response (data) scale
  fitted_response <- invlink(fitted_link, data_object$family, data_object$size)

  # analogous to the Gaussian case's SigInv_y - SigInv_X %*% betahat, but here
  # the latent link-scale vector w plays the role of y (GLMs work with a
  # Gaussian approximation on the link scale rather than on y directly)
  SigInv_r_list <- mapply(
    x = eigenprods_list, w = w_list, function(x, w) x$SigInv %*% as.matrix(w, ncol = 1) - x$SigInv_X %*% betahat,
    SIMPLIFY = FALSE
  )

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
    # same BLUP-averaging pattern as get_fitted_splm(): project onto each
    # level's indicator column, scale by the variance component, then average
    # nonzero contributions per level name
    fitted_randcov <- lapply(names(randcov_params), function(x) {
      fitted_val <- randcov_params[[x]] * do.call("rbind", mapply(
        z = data_object$randcov_list,
        r = SigInv_r_list,
        function(z, r) {
          crossprod(z[[x]][["Z"]], r)
        }
      ))
      fitted_val <- tapply(fitted_val, rownames(fitted_val), function(x) {
        if (any(x != 0)) {
          val <- mean(x[x != 0])
        } else {
          val <- 0
        }
      })
      # all combinations yields values with many zeros -- don't want to include these in the mean
      names_fitted_val <- rownames(fitted_val)
      fitted_val <- as.numeric(fitted_val)
      names(fitted_val) <- names_fitted_val
      fitted_val
    })
    names(fitted_randcov) <- names(randcov_params)
  }

  fitted_values <- list(
    response = fitted_response,
    link = fitted_link,
    spcov = list(de = fitted_de, ie = fitted_ie),
    randcov = fitted_randcov
  )
}

#' Get fitted values for an \code{spgautor()} model
#'
#' @param w The latent (link-scale) predictor vector
#' @param betahat Vector of fixed effects
#' @param spcov_params A \code{spcov_params} object
#' @param data_object The data object
#' @param eigenprods A \code{eigenprods} object
#' @param dist_matrix_list A list of distance matrices
#' @param randcov_params A \code{randcov_params} object
#'
#' @return A list of fitted values
#'
#' @noRd
get_fitted_spgautor <- function(w, betahat, spcov_params, data_object, eigenprods,
                                dist_matrix_list, randcov_params = NULL) {
  fitted_link <- as.numeric(w)
  # add offset
  if (!is.null(data_object$offset)) {
    fitted_link <- fitted_link + as.vector(data_object$offset)
  }
  fitted_response <- invlink(fitted_link, data_object$family, data_object$size)

  # subset the full neighborhood/weights structures down to just the observed
  # locations (spgautor models can include unobserved locations for prediction)
  dist_matrix <- data_object$W[data_object$observed_index, data_object$observed_index, drop = FALSE]
  M <- data_object$M[data_object$observed_index]

  # latent link-scale vector w plays the role of y in the Gaussian BLUP formula
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

  fitted_values <- list(
    response = fitted_response,
    link = fitted_link,
    spcov = list(de = fitted_de, ie = fitted_ie),
    randcov = fitted_randcov
  )
}

#' Apply the inverse link function for a GLM-type response family
#'
#' @param fitted_link Fitted values on the link scale
#' @param family The response family
#' @param size Binomial trial sizes (used only when \code{family} is \code{"binomial"})
#'
#' @return Fitted values on the response scale
#'
#' @noRd
invlink <- function(fitted_link, family, size) {
  # each branch is the inverse of that family's canonical link: exp() inverts
  # the log link (count/positive-continuous families), expit() (logistic
  # function) inverts the logit link (proportions), scaled by size for binomial counts
  if (family == "poisson") {
    fitted <- exp(fitted_link)
  } else if (family == "binomial") {
    if (is.null(size)) size <- 1
    fitted <- size * expit(fitted_link)
  } else if (family == "nbinomial") {
    fitted <- exp(fitted_link)
  } else if (family == "Gamma") {
    fitted <- exp(fitted_link)
  } else if (family == "inverse.gaussian") {
    fitted <- exp(fitted_link)
  } else if (family == "beta") {
    fitted <- expit(fitted_link)
  }
  fitted
}

#' Remove a model offset from a fitted link-scale vector
#'
#' @param w A fitted link-scale vector with the offset included, i.e.
#'   \code{fitted(object, type = "link")}
#' @param offset The model offset, or \code{NULL}
#'
#' @return \code{w} with the offset removed
#'
#' @details When a model has an offset, two distinct link-scale vectors are in
#'   play and using one where the other belongs produces plausible-looking but
#'   wrong numbers rather than an error. The offset-free latent vector
#'   \code{w = X beta + tau + epsilon} returned here is the process the spatial
#'   covariance describes, so it is the vector used by anything built from
#'   \code{Sigma}: kriging, the leave-one-out and k-fold updates, and the
#'   conditional-simulation residuals. The offset-inclusive linear predictor
#'   \code{w + offset} is the argument of the data model, so it is the vector
#'   used by anything built from the family: \code{get_d()}, \code{get_D()},
#'   \code{get_V()}, \code{get_var_y()}, \code{get_deviance_glm()}, and
#'   \code{invlink()}. Predictions are formed on the offset-free scale and the
#'   prediction location's own offset is added back at the end. 
#'
#' @noRd
w_offset_free <- function(w, offset) {
  if (is.null(offset)) {
    w
  } else {
    w - as.vector(offset)
  }
}

