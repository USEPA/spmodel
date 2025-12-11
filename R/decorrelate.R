#' Title
#'
#' @param formula filler
#' @param data filler
#' @param spcov_type filler
#' @param spcov_params filler
#' @param xcoord filler
#' @param ycoord filler
#' @param algorithm filler
#' @param training filler
#' @param evaluate_test filler
#' @param anisotropy filler
#' @param random filler
#' @param randcov_params filler
#' @param partition_factor filler
#' @param ordering filler
#' @param local filler
#' @param ... filler
#'
#' @return filler
#' @export
#' @references Matthew J. Heaton, Andrew Millane, and Jake S. Rhodes. 2025. A Scalable
#'   Spatial Decorrelation Preprocessing Approach for Machine and Deep Learning.
#'   \emph{Journal of Data Science}. 1-15, DOI 10.6339/25-JDS1210
#' @examples
#' out <- decorrelate(log_cond ~ temp, data = lake, spcov_type = "exponential")
#' out$grid
decorrelate <- function(formula, data, spcov_type, spcov_params, xcoord, ycoord, algorithm = "ranger", training, evaluate_test, anisotropy = FALSE, random, randcov_params, partition_factor, ordering = "maxmin", local, ...) {

  # set exponential as default if nothing specified
  if (missing(spcov_type) && missing(spcov_params)) {
    spcov_type <- "exponential"
    message("No spatial covariance type provided. Assuming \"exponential\".")
  }

  xcoord <- substitute(xcoord)
  ycoord <- substitute(ycoord)
  if (missing(spcov_params)) spcov_params <- NULL
  if (missing(random)) random <- NULL
  if (missing(randcov_params)) randcov_params <- NULL
  if (missing(partition_factor)) partition_factor <- NULL
  if (missing(local)) local <- NULL
  if (missing(training)) training <- NULL
  if (missing(evaluate_test)) evaluate_test <- FALSE


  if (is.null(spcov_params) || (!is.null(random) && is.null(randcov_params))) {
    evaluate_test <- TRUE
  }
  if (evaluate_test) {

    # if (missing(spcov_type)) spcov_type <- "exponential"

    init <- decorrelate_initial_search(
      formula = formula,
      data = data,
      spcov_type = spcov_type,
      spcov_params = spcov_params,
      xcoord = xcoord,
      ycoord = ycoord,
      algorithm = algorithm,
      training = training,
      anisotropy = anisotropy,
      random = random,
      randcov_params = randcov_params,
      partition_factor = partition_factor,
      ordering = ordering,
      local = local,
      ...
    )
    spcov_params <- init$spcov_params
    randcov_params <- init$randcov_params
    grid <- init$grid
    training_index <- init$training$training_index
    test_index <- init$training$test_index
    test_rmspe <- init$test_rmspe
  } else {
    grid <- NULL
    training_index <- NULL
    test_index <- NULL
    test_rmspe <- NULL
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
  obj <- list(
    algorithm = algorithm,
    decorrelate_data = decorr,
    fit = fit,
    grid = grid,
    test_index = test_index,
    training_index = training_index,
    test_rmspe = test_rmspe
  )
  new_obj <- structure(obj, class = "decorrelate")
  new_obj
}

fit_decorrelate_algorithm <- function(decorrelate_data, algorithm, ...){

  if (algorithm == "ranger") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Install the ranger package before using decorrelate() with algorithm \"ranger\".", call. = FALSE)
    }
    fit <- ranger::ranger(x = decorrelate_data$tX, y = decorrelate_data$ty, ...)
  } else if (algorithm == "randomForest") {
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("Install the randomForest package before using decorrelate() with algorithm \"randomForest\".", call. = FALSE)
    }
    fit <- randomForest::randomForest(x = decorrelate_data$tX, y = decorrelate_data$ty, ...)
  } else if (algorithm == "xgboost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("Install the xgboost package before using decorrelate() with algorithm \"xgboost\".", call. = FALSE)
    }
    fit <- xgboost::xgboost(x = decorrelate_data$tX, y = decorrelate_data$ty, ...)
  # } else if (algorithm == "nnet") {
  #   if (!requireNamespace("nnet", quietly = TRUE)) {
  #     stop("Install the nnet package before using decorrelate() with algorithm \"nnet\".", call. = FALSE)
  #   }
  #   fit <- nnet::nnet(x = decorrelate_data$tX, y = decorrelate_data$ty, ...)
  } else {
    stop("algorithm must be \"ranger\", \"randomForest\", or \"xgboost\". ", call. = FALSE)
  }

  fit

}

predict_decorrelate_algorithm <- function(decorrelate_data, tdata_test, algorithm, ...) {

  if (algorithm == "ranger") {
    preds <- predict(decorrelate_data, data = tdata_test$tX_newdata, ...)$predictions
  }

  if (algorithm == "randomForest") {
    preds <- predict(decorrelate_data, newdata = tdata_test$tX_newdata, ...)
  }

  if (algorithm == "xgboost") {
    preds <- predict(decorrelate_data, newdata = tdata_test$tX_newdata, ...)
  }

  if (algorithm == "nnet") {
    preds <- predict(decorrelate_data, newdata = tdata_test$tX_newdata, ...)
  }

  preds
}

get_training_list <- function(training, data) {

  if (is.null(training)) {
    training <- list(method = "split", prop = 0.75)
  }

  if (!is.list(training)) {
    stop("training must be a list")
  }

  names_training <- names(training)

  if (!"method" %in% names_training) {
    training$method <- "split"
  }

  if (!"prop" %in% names_training) {
    training$prop <- 0.75
  }

  if (!"training_index" %in% names_training && !"test_index" %in% names_training) {
    n <- NROW(data)
    index <- seq(1, n)
    index <- sample(index)
    n_training <- floor(n * training$prop)
    n_test <- n - n_training
    training$training_index <- index[seq(1, n_training)]
    training$test_index <- index[seq(n_training + 1, n)]
  } else if ("training_index" %in% names_training && !"test_index" %in% names_training) {
    n <- NROW(data)
    index <- seq(1, n)
    training$test_index <- index[-training$training_index]
  } else if (!"training_index" %in% names_training && "test_index" %in% names_training) {
    n <- NROW(data)
    index <- seq(1, n)
    training$training_index <- index[-training$test_index]
  }

  training$training_index <- sort(training$training_index)
  training$test_index <- sort(training$test_index)

  if (length(training$training_index) == 0) stop("No observations detected in training data.", call. = FALSE)
  if (length(training$test_index) == 0) stop("No observations detected in test data.", call. = FALSE)
  if (any(duplicated(training$training_index))) stop("Cannot have duplicated rows in training_index.", call. = FALSE)
  if (any(duplicated(training$test_index))) stop("Cannot have duplicated rows in test_index.", call. = FALSE)


  training

}
