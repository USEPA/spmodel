decorrelate_initial_search <- function(formula, data, spcov_type, spcov_params, xcoord, ycoord, algorithm, statistic, training, anisotropy, random, randcov_params, partition_factor, ordering = "maxmin", local, grid, ...) {


  training_list <- get_training_list(training, data)
  data_training <- data[training_list$training_index, , drop = FALSE]
  data_test <- data[training_list$test_index, , drop = FALSE]
  # yname <- as.character(attributes(terms(formula))$variables[[2]])
  yval <- model.response(model.frame(formula, data = data_test))


  grid_compare <- decorrelate_grid_internal(
    formula = formula,
    data = data_training,
    spcov_type = spcov_type,
    spcov_params = spcov_params,
    xcoord = xcoord,
    ycoord = ycoord,
    anisotropy = anisotropy,
    random = random,
    randcov_params = randcov_params
  )
  if (is.null(grid)) {
    grid <- grid_compare
  } else {
    if (any(! names(grid_compare) %in% names(grid))) stop("Invalid grid column names. Column names in grid must contain all column names returned by decorrelate_grid(formula, data, ...).", call. = FALSE)
    grid <- grid[, names(grid) %in% names(grid_compare), drop = FALSE]
  }

  params_list <- lapply(seq(1, NROW(grid)), function(x) {
    x <- grid[x, ]
    spcov_type <- x[["spcov_type"]]
    if ("extra" %in% names(x)) {
      spcov_params_val <- spcov_params(
        spcov_type = spcov_type,
        de = x[["de"]],
        ie = x[["ie"]],
        range = x[["range"]],
        extra = x[["extra"]],
        rotate = x[["rotate"]],
        scale = x[["scale"]]
      )
    } else {
      spcov_params_val <- spcov_params(
        spcov_type = spcov_type,
        de = x[["de"]],
        ie = x[["ie"]],
        range = x[["range"]],
        rotate = x[["rotate"]],
        scale = x[["scale"]]
      )
    }
    if (!is.null(random) || !is.null(randcov_params)) {
      remove_cols <- c("spcov_type", "de", "ie", "range", "rotate", "scale")
      if ("extra" %in% names(x)) remove_cols <- c(remove_cols, "extra")
      randcov_params_val <- unlist(x[, -which(names(x) %in% remove_cols), drop = FALSE])
      names(randcov_params_val) <- paste("(", names(randcov_params_val), ")", sep = "")
    } else {
      randcov_params_val <- NULL
    }
    list(spcov_params = spcov_params_val, randcov_params = randcov_params_val)
  })


  out <- lapply(params_list, function(x) {
    # warnings get repeated for each get_data_object() call
    tdata_training <- suppressWarnings(decorrelate_data_internal(
      formula = formula,
      data = data_training,
      spcov_params = x$spcov_params,
      xcoord = xcoord,
      ycoord = ycoord,
      randcov_params = x$randcov_params,
      partition_factor = partition_factor,
      ordering = ordering,
      local = local,
      ...
    ))
    # anisotropy is not getting accounted for somewhere here
    fit <- fit_decorrelate_algorithm(tdata_training, algorithm, ...)
    tdata_test <- decorrelate_newdata(tdata_training, newdata = data_test)
    preds <- predict_decorrelate_algorithm(fit, tdata_test, algorithm, ...)
    sp_decorr_preds <- recorrelate_newdata(tdata_test, preds)
    errors <- yval - sp_decorr_preds
    rmspe <- sqrt(mean(errors^2))
    medae <- median(abs(errors))
    cor2 <- cor(yval, sp_decorr_preds)^2
    list(rmspe = rmspe, medae = medae, cor2 = cor2)
  })
  grid$rmspe <- unlist(lapply(out, function(x) x$rmspe))
  grid$medae <- unlist(lapply(out, function(x) x$medae))
  grid$cor2 <- unlist(lapply(out, function(x) x$cor2))
  if (statistic == "cor2") {
    best_val <- which.max(grid[[statistic]])
  } else {
    best_val <- which.min(grid[[statistic]])
  }

  test_rmspe <- grid$rmspe[best_val]
  test_medae <- grid$medae[best_val]
  test_cor2 <- grid$cor2[best_val]
  spcov_params_val <- params_list[[best_val]]$spcov_params
  randcov_params_val <- params_list[[best_val]]$randcov_params

  if (statistic == "cor2") {
    grid <- grid[order(grid[[statistic]], decreasing = TRUE), , drop = FALSE]
  } else {
    grid <- grid[order(grid[[statistic]]), , drop = FALSE]
  }
  row.names(grid) <- NULL

  list(spcov_params = spcov_params_val, randcov_params = randcov_params_val, grid = grid, training = training_list, test_rmspe = test_rmspe, test_medae = test_medae, test_cor2 = test_cor2)
}
