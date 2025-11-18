decorrelate_newdata <- function(object, newdata, local, ...) {


  if (!inherits(object, "decorrelate")) {
    stop("object must have class \"decorrelate\".", call. = FALSE)
  }

  if (!missing(local)) {
    object$local <- get_local_list_decorrelate(local)
  }

  # rename relevant quantities
  obdata <- object$obdata
  xcoord <- object$xcoord
  ycoord <- object$ycoord

  if (missing(newdata)) {
    newdata <- object$newdata
  }

  # error if newdata missing from arguments and object
  if (missing(newdata) && is.null(object$newdata)) {
    stop("No missing data to predict. newdata must be specified in the newdata argument or object$newdata must be non-NULL.", call. = FALSE)
  }

  # save spcov param vector
  spcov_params_val <- object$coefficients$spcov

  # save randcov param vector
  randcov_params_val <- object$coefficients$randcov

  # partition factor
  partition_factor_val <- object$partition_factor

  attr_sp <- attr(class(newdata), "package")
  if (!is.null(attr_sp) && length(attr_sp) == 1 && attr_sp == "sp") {
    stop("sf objects must be used instead of sp objects. To convert your sp object into an sf object, run sf::st_as_sf().", call. = FALSE)
  }

  if (inherits(newdata, "sf")) {
    newdata <- suppressWarnings(sf::st_centroid(newdata))

    newdata <- sf_to_df(newdata)
    names(newdata)[[which(names(newdata) == ".xcoord")]] <- as.character(xcoord) # only relevant if newdata is sf data is not
    names(newdata)[[which(names(newdata) == ".ycoord")]] <- as.character(ycoord) # only relevant if newdata is sf data is not
  }

  # add back in zero column to cover anisotropy (should make anisotropy only available 1-d)
  if (object$dim_coords == 1) {
    obdata[[ycoord]] <- 0
    newdata[[ycoord]] <- 0
  }

  if (object$anisotropy) { # could just do rotate != 0 || scale != 1
    obdata_aniscoords <- transform_anis(obdata, xcoord, ycoord,
                                        rotate = spcov_params_val[["rotate"]],
                                        scale = spcov_params_val[["scale"]]
    )
    obdata[[xcoord]] <- obdata_aniscoords$xcoord_val
    obdata[[ycoord]] <- obdata_aniscoords$ycoord_val
    object$obdata <- obdata
    newdata_aniscoords <- transform_anis(newdata, xcoord, ycoord,
                                         rotate = spcov_params_val[["rotate"]],
                                         scale = spcov_params_val[["scale"]]
    )
    newdata[[xcoord]] <- newdata_aniscoords$xcoord_val
    newdata[[ycoord]] <- newdata_aniscoords$ycoord_val
  }

  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  if (any(grepl("nmatrix.", attributes(formula_newdata)$dataClasses, fixed = TRUE)) && NROW(newdata) == 1) {
    newdata <- newdata[c(1, 1), , drop = FALSE]
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
    newdata_model <- newdata_model[1, , drop = FALSE]
    # find offset
    offset <- model.offset(newdata_model_frame)
    if (!is.null(offset)) {
      offset <- offset[1]
    }
    newdata <- newdata[1, , drop = FALSE]
  } else {
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    # assumes that predicted observations are not outside the factor levels
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
  }
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(object$X)) # colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # finding rows w/out NA
  ob_predictors <- complete.cases(newdata_model)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }

  # storing newdata as a list
  newdata_rows_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_list <- mapply(x = newdata_rows_list, y = newdata_model_list, FUN = function(x, y) list(row = x, x0 = y), SIMPLIFY = FALSE)

  # randcov stuff
  extra_randcov_list <- get_extra_randcov_list(object, obdata, newdata)
  # reform_bar2_list <- extra_randcov_list$reform_bar2_list
  # Z_index_obdata_list <- extra_randcov_list$Z_index_obdata_list
  # reform_bar1_list <- extra_randcov_list$reform_bar1_list
  # Z_val_obdata_list <- extra_randcov_list$Z_val_obdata_list

  # partition stuff
  extra_partition_list <- get_extra_partition_list(object, obdata, newdata)
  # reform_bar2 <- extra_partition_list$reform_bar2
  # partition_index_obdata <- extra_partition_list$partition_index_obdata

  if (object$local$method == "all") {
    cov_mat <- covmatrix.splm(object)
    cor_mat <- cov_mat / object$total_var
    cor_lowchol <- t(chol(cor_mat))
    rSqrtSigInv_X <- forwardsolve(cor_lowchol, object$X)
    rSqrtSigInv_y <- forwardsolve(cor_lowchol, object$y)
    cor_lowchol_list <- list(
      cor_lowchol = cor_lowchol,
      rSqrtSigInv_X = rSqrtSigInv_X,
      rSqrtSigInv_y = rSqrtSigInv_y
    )
  } else {
    cor_lowchol_list <- NULL
  }

  if (object$local$parallel) {
    cl <- parallel::makeCluster(object$local$ncores)
    output <- parallel::parLapply(cl, newdata_list, get_decorrelate_newdata,
                                  object, cor_lowchol_list, extra_randcov_list, extra_partition_list)
    cl <- parallel::stopCluster(cl)
  } else {
    output <- lapply(newdata_list, get_decorrelate_newdata,
                                  object, cor_lowchol_list, extra_randcov_list, extra_partition_list)
  }

  tX_newdata <- do.call("rbind", lapply(output, function(x) x$tX_newdata))
  rownames(tX_newdata) <- rownames(newdata)
  colnames(tX_newdata) <- colnames(object$X)
  yscale <- do.call("c", lapply(output, function(x) x$yscale))
  names(yscale) <- rownames(newdata)
  yoffset <- do.call("c", lapply(output, function(x) x$yoffset))
  names(yoffset) <- rownames(newdata)
  output <- list(
    tX_newdata = tX_newdata,
    yscale = yscale,
    yoffset = yoffset
  )
  new_output <- structure(output, class = "decorrelate_newdata")
  new_output

}

get_decorrelate_newdata <- function(newdata_list, object, cor_lowchol_list, extra_randcov_list, extra_partition_list) {

  obdata <- object$obdata
  xcoord <- object$xcoord
  ycoord <- object$ycoord
  X <- object$X
  y <- object$y

  # storing partition vector
  if (!is.null(object$partition_factor)) {
    partition_vector <- partition_vector(object$partition_factor,
                                         data = object$obdata,
                                         newdata = newdata_list$row, reform_bar2 = extra_partition_list$reform_bar2,
                                         partition_index_data = extra_partition_list$partition_index_obdata
    )
  } else {
    partition_vector <- NULL
  }

  dist_vector <- spdist_vectors(newdata_list$row, obdata, xcoord, ycoord, object$dim_coords)

  # making random vector if necessary
  if (!is.null(object$random)) {
    randcov_vector_val <- randcov_vector(object$coefficients$randcov, object$obdata, newdata_list$row,
                                         extra_randcov_list$reform_bar2_list, extra_randcov_list$Z_index_obdata_list)
  } else {
    randcov_vector_val <- NULL
  }

  # making the covariance vector
  cov_vector_val <- cov_vector(object$coefficients$spcov, dist_vector, randcov_vector_val, partition_vector)
  cov_vector_val <- as.numeric(cov_vector_val)

  # subsetting data if method distance
  if (object$local$method == "distance") {
    n <- length(cov_vector_val)
    # want the smallest distance here and order goes from smallest first to largest last (keep last values with are smallest distance)
    nn_index <- order(as.numeric(dist_vector))[seq(from = 1, to = min(n, object$local$size))]
    obdata <- obdata[nn_index, , drop = FALSE]
    X <- X[nn_index, , drop = FALSE]
    y <- y[nn_index]
    cov_vector_val <- cov_vector_val[nn_index]
  }

  if (object$local$method == "covariance") {
    n <- length(cov_vector_val)
    # want the largest covariance here and order goes from smallest first to largest last (keep last values which are largest covariance)
    cov_index <- order(as.numeric(cov_vector_val))[seq(from = n, to = max(1, n - object$local$size + 1))] # use abs() here?
    obdata <- obdata[cov_index, , drop = FALSE]
    X <- X[cov_index, , drop = FALSE]
    y <- y[cov_index]
    cov_vector_val <- cov_vector_val[cov_index]
  }

  if (object$local$method %in% c("distance", "covariance")) {
    if (!is.null(object$random)) {
      randcov_names <- get_randcov_names(object$random)
      xlev_list <- lapply(extra_randcov_list$Z_index_obdata_list, function(x) x$reform_bar2_xlev)
      randcov_Zs <- get_randcov_Zs(obdata, randcov_names, xlev_list = xlev_list)
    }
    partition_matrix_val <- partition_matrix(object$partition_factor, obdata)
    cov_matrix_val <- cov_matrix(
      object$coefficients$spcov, spdist(obdata, xcoord, ycoord), object$coefficients$randcov,
      randcov_Zs, partition_matrix_val,
      diagtol = object$diagtol
    )
    cor_matrix_val <- cov_matrix_val / object$total_var
    cor_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(cor_matrix_val)))
  } else {
    cor_lowchol <- cor_lowchol_list$cor_lowchol
  }

  cor_vector_val <- cov_vector_val / object$total_var

  rSqrtSigInv_r0 <- forwardsolve(cor_lowchol, cor_vector_val)
  r0_SigInv_r0 <- crossprod(rSqrtSigInv_r0, rSqrtSigInv_r0)
  w <- as.numeric(1 - r0_SigInv_r0)

  if (object$local$method %in% c("distance", "covariance")) {
    rSqrtSigInv_X <- forwardsolve(cor_lowchol, X)
    rSqrtSigInv_y <- forwardsolve(cor_lowchol, y)
  } else {
    rSqrtSigInv_X <- cor_lowchol_list$rSqrtSigInv_X
    rSqrtSigInv_y <- cor_lowchol_list$rSqrtSigInv_y
  }

  sqrt_w <- sqrt(w)
  tX_newdata <- (newdata_list$x0 - crossprod(rSqrtSigInv_r0, rSqrtSigInv_X)) / sqrt_w
  yoffset <- crossprod(rSqrtSigInv_r0, rSqrtSigInv_y)

  list(
    tX_newdata = tX_newdata,
    yscale = sqrt_w,
    yoffset = as.numeric(yoffset)
  )


}
