decorrelate <- function(formula, data, spcov_type, spcov_params, xcoord, ycoord, algorithm = "ranger", train, anisotropy = FALSE, randcov_params, partition_factor, ordering = "maxmin", local, ...) {

  xcoord <- substitute(xcoord)
  ycoord <- substitute(ycoord)
  if (missing(randcov_params)) randcov_params <- NULL
  if (missing(partition_factor)) partition_factor <- NULL
  if (missing(local)) local <- NULL
  if (missing(train)) train <- NULL

  if (missing(spcov_params)) {

    if (missing(spcov_type)) spcov_type <- "exponential"
    init <- decorrelate_initial_search(
      formula = formula,
      data = data,
      spcov_type = spcov_type,
      xcoord = xcoord,
      ycoord = ycoord,
      algorithm = algorithm,
      train = train,
      anisotropy = anisotropy,
      randcov_params = randcov_params,
      partition_factor = partition_factor,
      ordering = ordering,
      local = local,
      ...
    )
    spcov_params <- init$spcov_params
    randcov_params <- init$randcov_params
  }

  decorr <- decorrelate_data(
    formula = formula,
    data = data,
    spcov_params = spcov_params,
    xcoord = xcoord,
    ycoord = ycoord,
    randcov_params = randcov_params,
    partition_factor = partition_factor,
    ordering = ordering,
    local = local,
    ...
  )

  fit <- fit_decorrelate_algorithm(decorr, algorithm, ...)
  obj <- list(decorrelate_data = decorr, fit = fit, algorithm = algorithm, test_rmspe = init$min_rmspe)
  new_obj <- structure(obj, class = "decorrelate")
  new_obj
}

fit_decorrelate_algorithm <- function(object, algorithm, ...){

  if (algorithm == "ranger") {

    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Install the ranger package before using decorrelate() with algorithm \"ranger\".", call. = FALSE)
    }

    fit <- ranger::ranger(x = object$tX, y = object$ty, ...)

  }

  # xgboost, nnet, randomForest

  fit

}

predict_decorrelate_algorithm <- function(object, tdata_test, algorithm, ...) {

  if (algorithm == "ranger") {

    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Install the ranger package before using decorrelate() with algorithm \"ranger\".", call. = FALSE)
    }

    preds <- predict(object, data = tdata_test$tX_newdata, ...)$predictions

  }

  # xgboost, nnet, randomForest

  preds
}

get_train_list <- function(train, data) {

  if (is.null(train)) {
    train <- list(method = "split", prop = 0.75)
  }

  if (!is.list(train)) {
    stop("train must be a list")
  }

  names_train <- names(train)

  if (!"method" %in% names_train) {
    train$method <- "split"
  }

  if (!"prop" %in% names_train) {
    train$prop <- 0.75
  }

  if (!"train_index" %in% names_train) {
    n <- NROW(data)
    n_train <- floor(n * train$prop)
    n_test <- n - n_train
    index <- sample(seq(1, n))
    train$train_index <- sort(index[seq(1, n_train)])
    train$test_index <- sort(index[seq(n_train + 1, n)])
  }

  if (!"test_index" %in% names_train) {
    n <- NROW(data)
    n_train <- floor(n * train$prop)
    n_test <- n - n_train
    index <- seq(1, n)
    train$train_index <- sort(train$train_index)
    train$test_index <- sort(index[-train$train_index])
  }


  train

}
