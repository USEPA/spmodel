get_cov_gradients_context_splm <- function(object) {

  spcov_params_val <- coef(object, type = "spcov")
  spcov_type <- class(spcov_params_val)
  spcov_is_known <- object$is_known$spcov

  spcov_initial_val <- do.call(
    spcov_initial,
    c(
      list(spcov_type = spcov_type),
      as.list(spcov_params_val),
      list(known = names(spcov_is_known)[spcov_is_known])
    )
  )

  has_randcov <- !is.null(object$random)
  if (has_randcov) {
    randcov_params_val <- coef(object, type = "randcov")
    randcov_is_known <- object$is_known$randcov
    randcov_initial_val <- do.call(
      randcov_initial,
      c(
        as.list(randcov_params_val),
        list(known = names(randcov_is_known)[randcov_is_known])
        )
      )
  } else {
    randcov_params_val <- NULL
    randcov_is_known <- NULL
    randcov_initial_val <- NULL
  }

  data_object <- get_data_object_splm(
    formula = object$formula, data = object$obdata, spcov_initial = spcov_initial_val,
    xcoord = object$xcoord, ycoord = object$ycoord, estmethod = object$estmethod,
    anisotropy = object$anisotropy, random = object$random, randcov_initial = randcov_initial_val,
    partition_factor = object$partition_factor, local = FALSE, range_constrain = FALSE
  )

  X <- data_object$X_list[[1]]
  diagtol <- object$diagtol
  randcov_Zs <- if (has_randcov) data_object$randcov_list[[1]] else NULL
  partition_matrix_val <- if (!is.null(data_object$partition_list)) data_object$partition_list[[1]] else NULL

  if (object$anisotropy) {
    dist_matrix_list <- NULL
    dist_matrix <- NULL
  } else if (inherits(spcov_params_val, c("none", "ie")) && !has_randcov) {
    dist_matrix_list <- NULL
    dist_matrix <- NULL
  } else {
    dist_matrix_list <- lapply(data_object$obdata_list, function(d) spdist(d, data_object$xcoord, data_object$ycoord))
    dist_matrix <- as.matrix(dist_matrix_list[[1]])
  }


  # deprofiling step
  spcov_orig2optim_val <- spcov_orig2optim(spcov_initial_val, spcov_profiled = FALSE, data_object = data_object)
  randcov_orig2optim_val <- if (has_randcov) {
    randcov_orig2optim(randcov_initial_val, randcov_profiled = FALSE, spcov_initial = spcov_initial_val)
  } else {
    NULL
  }
  eta_val <- assemble_optim_par(spcov_orig2optim_val, randcov_orig2optim_val)

  spcov_names_free <- names(spcov_is_known)[!spcov_is_known]
  randcov_names_free <- if (has_randcov) names(randcov_is_known)[!randcov_is_known] else character(0)
  cov_names_free <- c(spcov_names_free, randcov_names_free)

  if (length(cov_names_free) == 0) {
    stop("All covariance parameters known. Satterthwaite not applicable.", call. = FALSE)
  }

  cov_val <- c(as.numeric(spcov_params_val), as.numeric(randcov_params_val))
  names(cov_val) <- c(names(spcov_params_val), names(randcov_params_val))
  cov_val_free <- cov_val[cov_names_free]

  context <- list(
    data_object = data_object, X = X, diagtol = diagtol,
    anisotropy = object$anisotropy, estmethod = object$estmethod,
    spcov_params = spcov_params_val, spcov_is_known = spcov_is_known,
    randcov_params = randcov_params_val, randcov_is_known = randcov_is_known, randcov_Zs = randcov_Zs,
    spcov_names_free = spcov_names_free,
    randcov_names_free = randcov_names_free, cov_names_free = cov_names_free,
    cov_val_free = cov_val_free,
    partition_matrix = partition_matrix_val,
    dist_matrix = dist_matrix,
    spcov_orig2optim = spcov_orig2optim_val, randcov_orig2optim = randcov_orig2optim_val,
    eta = eta_val
  )


}
