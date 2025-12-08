decorrelate_initial_search <- function(formula, data, spcov_type, xcoord, ycoord, algorithm, train, anisotropy, randcov_params, partition_factor, ordering = "maxmin", local, ...) {

  if (!is.null(randcov_params)) {
    random <- reformulate(names(randcov_params))
  } else {
    random <- NULL
  }


  grid <- decorrelate_grid(
    formula = formula,
    data = data,
    spcov_type = spcov_type,
    xcoord = xcoord,
    ycoord = ycoord,
    anisotropy = anisotropy,
    random = random
  )

  if (!anisotropy) {
    grid$rotate <- 0
    grid$scale <- 1
  }

  params_list <- lapply(seq(1, NROW(grid)), function(x) {
    x <- grid[x, ]
    spcov_params_val <- spcov_params(
      spcov_type = x[["spcov_type"]],
      de = x[["de"]],
      ie = x[["ie"]],
      range = x[["range"]],
      rotate = x[["rotate"]],
      scale = x[["scale"]]
    )
    if (!is.null(random)) {
      remove_cols <- c("spcov_type", "de", "ie", "range", "rotate", "scale")
      randcov_params_val <- unlist(x[, -which(names(x) %in% remove_cols), drop = FALSE])
      names(randcov_params_val) <- paste("(", names(randcov_params_val), ")", sep = "")
    } else {
      randcov_params_val <- NULL
    }
    list(spcov_params = spcov_params_val, randcov_params = randcov_params_val)
  })

  train_list <- get_train_list(train, data)
  data_train <- data[train_list$train_index, , drop = FALSE]
  data_test <- data[train_list$test_index, , drop = FALSE]
  yname <- as.character(attributes(terms(formula))$variables[[2]])
  # need to come back and specify x levels in the training data here

  out <- lapply(params_list, function(x) {
    tdata_train <- decorrelate_data(
      formula = formula,
      data = data_train,
      spcov_params = x$spcov_params,
      randcov_params = x$randcov_params,
      partition_factor = partition_factor,
      ordering = ordering,
      local = local,
      ...
    )
    fit <- fit_decorrelate_algorithm(tdata_train, algorithm, ...)
    tdata_test <- decorrelate_newdata(tdata_train, newdata = data_test)
    preds <- predict_decorrelate_algorithm(fit, tdata_test, algorithm, ...)
    sp_decorr_preds <- recorrelate_newdata(preds, tdata_test)
    errors <- data_test[[yname]] - sp_decorr_preds
    rmspe <- sqrt(mean(errors^2))
    rmspe
  })
  min_rmspe <- which.min(unlist(out))
  list(spcov_params = params_list[[min_rmspe]]$spcov_params, randcov_params = params_list[[min_rmspe]]$randcov_params, min_rmspe = out[[min_rmspe]])
}
