#' Apply the Spatial Decorrelation Transformation for Machine Learning Models
#'
#' @description Apply the spatial decorrelation transformation for
#'   point-referenced data, allowing for random effects,
#'   anisotropy, partition factors, and big data methods.
#'
#' @param formula A two-sided linear formula describing the fixed effect structure
#'   of the model, with the response to the left of the \code{~} operator and
#'   the terms on the right, separated by \code{+} operators.
#' @param data A data frame or \code{sf} object object that contains
#'   the variables in \code{fixed}, \code{random}, and \code{partition_factor}
#'   as well as geographical information. If an \code{sf} object is
#'   provided with \code{POINT} geometries, the x-coordinates and y-coordinates
#'   are used directly. If an \code{sf} object is
#'   provided with \code{POLYGON} geometries, the x-coordinates and y-coordinates
#'   are taken as the centroids of each polygon.
#' @param spcov_type The spatial covariance type. Available options include
#'   \code{"exponential"}, \code{"spherical"}, \code{"gaussian"},
#'   \code{"triangular"}, \code{"circular"}, \code{"cubic"},
#'   \code{"pentaspherical"}, \code{"cosine"}, \code{"wave"},
#'   \code{"jbessel"}, \code{"gravity"}, \code{"rquad"},
#'   \code{"magnetic"}, \code{"matern"}, \code{"cauchy"}, \code{"pexponential"},
#'   and \code{"none"}. Parameterizations of each spatial covariance type are
#'   available in Details. Multiple spatial covariance types can be provided as
#'   a character vector, and then \code{decorrelate()} is called iteratively for each
#'   element and a list is returned for each model fit. The default for
#'   \code{spcov_type} is \code{"exponential"}. When \code{spcov_type} is
#'   specified, all spatial covariance parameters are estimated.
#'   \code{spcov_type} is ignored if \code{spcov_params} is provided.
#' @param spcov_params An object from [spcov_params()] that contains the
#'   spatial covariance parameters used by the spatial decorrelation transformation.
#' @param xcoord The name of the column in \code{data} representing the x-coordinate.
#'   Can be quoted or unquoted. Not required if \code{data} is an \code{sf} object.
#' @param ycoord The name of the column in \code{data} representing the y-coordinate.
#'   Can be quoted or unquoted. Not required if \code{data} is an \code{sf} object.
#' @param algorithm The machine learning algorithm applied. Available options
#'   include \code{"ranger"}, \code{"randomForest"}, and \code{"xgboost"}.
#'   \code{"ranger"} specifies a random forest via [ranger::ranger()].
#'   \code{"randomForest"} specifies a random forest via [randomForest::randomForest()].
#'   \code{"xgboost"} specifies a boosted decision tree ensemble via [xgboost::xgboost()].
#' @param training An list controlling how the training and test data are assigned
#'   when evaluating test data performance.
#'   The following arguments detail this process:
#'   \itemize{
#'    \item \code{method: }{The method used to evaluate test data performance.
#'      Currently, the only option is \code{"split"}, which splits \code{data} up
#'      into distinct training and test sets.}
#'    \item \code{prop: }{The proportion (a numeric vector between zero and one) of observations in \code{data} that should
#'      be assigned to the training data. The default is 0.75, which means that
#'      75\% of the observations are assigned to the training data and 25% to the
#'      test data. Ignored if \code{training_index} or \code{test_index} are provided.}
#'    \item \code{training_index: }{A numeric vector that specifies which rows (i.e., indices)
#'      of \code{data} should be assigned to the training data. If omitted, defaults
#'      to the rows which are not already included in \code{test_index}.}
#'    \item \code{test_index: }{A numeric vector that specifies which rows (i.e., indices)
#'      of \code{data} should be assigned to the test data. If omitted, defaults
#'      to the rows which are not already included in \code{training_index}.}
#'   }
#'   If omitted, \code{training} is transformed into
#'   \code{list(method = "split", prop = 0.75)}.
#' @param evaluate_test A logical indicating whether a grid should be constructed
#'   and evaluated when spatial decorrelation parameters are known (i.e.,
#'   \code{spcov_params} is specified, and, if random effects are included, \code{randcov_params} is specified).
#'   If \code{TRUE}, constructs and evalutes the grid after assigning observations to
#'   test and training data sets. If any parameters (spatial or random effects) are estimated via a grid search,
#'   \code{evaluate_test} is set to \code{TRUE}.
#' @param anisotropy A logical indicating whether (geometric) anisotropy should
#'   be modeled. Not required if the \code{rotate} and \code{scale} parameters in \code{spcov_params()} are
#'   0 and 1, respectively. When \code{anisotropy} is \code{TRUE},
#'   computational times can significantly increase. The default is \code{FALSE}.
#' @param random A one-sided linear formula describing the random effect structure
#'   of the model. Terms are specified to the right of the \code{~ operator}.
#'   Each term has the structure \code{x1 + ... + xn | g1/.../gm}, where \code{x1 + ... + xn}
#'   specifies the model for the random effects and \code{g1/.../gm} is the grouping
#'   structure. Separate terms are separated by \code{+} and must generally
#'   be wrapped in parentheses. Random intercepts are added to each model
#'   implicitly when at least  one other variable is defined.
#'   If a random intercept is not desired, this must be explicitly
#'   defined (e.g., \code{x1 + ... + xn - 1 | g1/.../gm}). If only a random intercept
#'   is desired for a grouping structure, the random intercept must be specified
#'   as \code{1 | g1/.../gm}. Note that \code{g1/.../gm} is shorthand for \code{(1 | g1/.../gm)}.
#'   If only random intercepts are desired and the shorthand notation is used,
#'   parentheses can be omitted.
#' @param randcov_params An object from [randcov_params()] that contains the
#'   random effect variances used by the spatial decorrelation transformation.
#' @param partition_factor A one-sided linear formula with a single term
#'   specifying the partition factor.  The partition factor assumes observations
#'   from different levels of the partition factor are uncorrelated.
#' @param ordering The data ordering applied. Available options
#'   include \code{"maxmin"}, \code{"grts"}, \code{"random"}, and \code{"xgboost"}.
#'   \code{"maxmin"} applies maximum minimum distance ordering via \code{GPvecchia::order_maxmin_exact()}.
#'   \code{"maxmin"} applies ordering using a spatially balanced GRTS sample via \code{spsurvey::grts()}.
#'   \code{"random"} applies a completely random ordering.
#'   \code{"maxmin"} applies no random ordering.
#'   The default is \code{"maxmin"}.
#' @param local A optional logical or list controlling the big data approximation.
#'   If omitted, \code{local} is set
#'   to \code{TRUE} or \code{FALSE} based on the sample size (the number of
#'   non-missing observations in \code{data}) -- if the sample size exceeds 5,000,
#'   \code{local} is set to \code{TRUE}. Otherwise it is set to \code{FALSE}.
#'   If \code{local} is \code{FALSE}, no big data approximation
#'   is implemented. If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{method}: The big data approximation method. If \code{method = "all"},
#'       all observations are used and \code{size} is ignored. If \code{method = "distance"},
#'       the \code{size} data observations closest (in terms of Euclidean distance)
#'       to the observation requiring prediction are used.
#'       If \code{method = "covariance"}, the \code{size} data observations
#'       with the highest covariance with the observation requiring prediction are used.
#'       If random effects and partition factors are not used in estimation and
#'       the spatial covariance function is monotone decreasing,
#'       \code{"distance"} and \code{"covariance"} are equivalent. The default
#'       is \code{"covariance"}.
#'     \item \code{size}: The number of data observations to use when \code{method}
#'       is \code{"distance"} or \code{"covariance"}. The default is 30.
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. This can significantly speed
#'       up computations even when \code{method = "all"} (i.e., no big data
#'       approximation is used), as predictions
#'       are spread out over multiple cores. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 30, method = "covariance", parallel = FALSE)}.
#' @param ... Other arguments to the functions called by \code{algorithm}.
#'
#' @details The spatial decorrelation transformation is a preprocessing transformation
#'   that reduces the impacts of spatial dependence (i.e., covariance, correlation)
#'   on machine learning models. Predictions are made on the
#'   decorrelated scale and then recorrelated to account for spatial dependence.
#'   See Heaton et al., 2025 for details.
#'
#'   \code{spcov_type} Details: The correlation matrix \eqn{R} controls the spatial dependence structure
#'   among observations. Parametric forms for \eqn{R} are given below, where \eqn{\eta = h / range}
#'   for \eqn{h} distance between observations:
#'   \itemize{
#'     \item exponential: \eqn{exp(- \eta )}
#'     \item spherical: \eqn{(1 - 1.5\eta + 0.5\eta^3) * I(h <= range)}
#'     \item gaussian: \eqn{exp(- \eta^2 )}
#'     \item triangular: \eqn{(1 - \eta) * I(h <= range)}
#'     \item circular: \eqn{(1 - (2 / \pi) * (m * sqrt(1 - m^2) + sin^{-1}(m))) * I(h <= range), m = min(\eta, 1)}
#'     \item cubic: \eqn{(1 - 7\eta^2 + 8.75\eta^3 - 3.5\eta^5 + 0.75\eta^7) * I(h <= range)}
#'     \item pentaspherical: \eqn{(1 - 1.875\eta + 1.25\eta^3 - 0.375\eta^5) * I(h <= range)}
#'     \item cosine: \eqn{cos(\eta)}
#'     \item wave: \eqn{sin(\eta) / \eta * I(h > 0) + I(h = 0)}
#'     \item jbessel: \eqn{Bj(h * range)}, Bj is Bessel-J function
#'     \item gravity: \eqn{(1 + \eta^2)^{-0.5}}
#'     \item rquad: \eqn{(1 + \eta^2)^{-1}}
#'     \item magnetic: \eqn{(1 + \eta^2)^{-1.5}}
#'     \item matern: \eqn{2^{1 - extra}/ \Gamma(extra) * \alpha^{extra} * Bk(\alpha, extra)}, \eqn{\alpha = (2extra * \eta)^{0.5}}, Bk is Bessel-K function with order \eqn{1/5 \le extra \le 5}
#'     \item cauchy: \eqn{(1 + \eta^2)^{-extra}}, \eqn{extra > 0}
#'     \item pexponential: \eqn{exp(h^{extra}/range)}, \eqn{0 < extra \le 2}
#'     \item none: \eqn{0}
#'   }
#'
#'   All spatial covariance functions are valid in one spatial dimension. All
#'   spatial covariance functions except \code{triangular} and \code{cosine} are
#'   valid in two dimensions. An alias for \code{none} is \code{ie}.
#'
#' \code{anisotropy} Details: By default, all spatial covariance parameters except \code{rotate}
#'   and \code{scale} as well as all random effect variance parameters
#'   are assumed unknown, requiring estimation. If either \code{rotate} or \code{scale}
#'   are given initial values other than 0 and 1 (respectively)
#'   in [spcov_params()], \code{anisotropy} is implicitly set to \code{TRUE}.
#'   (Geometric) Anisotropy is modeled by transforming a covariance function that
#'   decays differently in different directions to one that decays equally in all
#'   directions via rotation and scaling of the original coordinates. The rotation is
#'   controlled by the \code{rotate} parameter in \eqn{[0, \pi]} radians. The scaling
#'   is controlled by the \code{scale} parameter in \eqn{[0, 1]}. The anisotropy
#'   correction involves first a rotation of the coordinates clockwise by \code{rotate} and then a
#'   scaling of the coordinates' minor axis by the reciprocal of \code{scale}. The spatial
#'   covariance is then computed using these transformed coordinates.
#'
#'  \code{random} Details: If random effects are used (the estimation method must be \code{"reml"} or
#'   \code{"ml"}), the model
#'   can be written as \eqn{y = X \beta + Z1u1 + ... Zjuj + \tau + \epsilon},
#'   where each Z is a random effects design matrix and each u is a random effect.
#'
#'  \code{partition_factor} Details: The partition factor can be represented in matrix form as \eqn{P}, where
#'   elements of \eqn{P} equal one for observations in the same level of the partition
#'   factor and zero otherwise. The covariance matrix involving only the
#'   spatial and random effects components is then multiplied element-wise
#'   (Hadmard product) by \eqn{P}, yielding the final covariance matrix.
#'
#'   \code{local} Details: The big data approximation works by leveraging the
#'   conditional nature of the spatial decorrelation transformation via the
#'   Vecchia approximation. The Vecchia approximation enables efficient computation
#'   of the conditional distribution by considering only the \code{size} most relevant
#'   observations in the ordering (rather than using all the observations).
#'
#'   Observations with \code{NA} response values are removed for model
#'   fitting, but their values can be predicted afterwards by running
#'   \code{predict(object)}.
#'
#' @return A list with many elements that store information about the fitted model object:
#'   \itemize{
#'     \item{algorithm: }{The machine learning algorithm used.}
#'     \item{decorrelate_data: }{The output of [decorrelate_data()] applied to \code{data}.}
#'     \item{fit: }{The fitted machine learning model object applied to the decorrelated data.}
#'     \item{grid: }{If used, the grid of spatial decorrelation parameters evaluated and their corresponding
#'       root-mean-squared-prediction error (rmspe) when applied to the test data.}
#'     \item{training_index: }{If used, the rows (i.e., indices) in \code{data} that were assigned to the training data.}
#'     \item{test_index: }{If used, the rows (i.e., indices) in \code{data} that were assigned to the test data.}
#'   }
#'
#' @export
#'
#' @references Matthew J. Heaton, Andrew Millane, and Jake S. Rhodes. 2025. A Scalable
#'   Spatial Decorrelation Preprocessing Approach for Machine and Deep Learning.
#'   \emph{Journal of Data Science}. 1-15, DOI 10.6339/25-JDS1210
#'
#' @examples
#' out <- decorrelate(log_cond ~ temp, data = lake, spcov_type = "exponential")
#' out$grid
decorrelate <- function(formula, data, spcov_type, spcov_params, xcoord, ycoord, algorithm = "ranger", training, evaluate_test, anisotropy = FALSE, random, randcov_params, partition_factor, ordering = "maxmin", local, ...) {

  # set exponential as default if nothing specified
  if (missing(spcov_type) && missing(spcov_params)) {
    spcov_type <- "exponential"
    message("No spatial covariance type provided. Assuming \"exponential\".")
  }

  if (!missing(spcov_type) && length(spcov_type) > 1) {
    if (missing(training)) training <- list()
    call_list <- as.list(match.call())[-1]
    call_list$training <- get_training_list(training, data)
    penv <- parent.frame()
    decorrelate_out <- lapply(spcov_type, function(x) {
      call_list$spcov_type <- x
      do.call("decorrelate", call_list, envir = penv)
    })
    names(decorrelate_out) <- spcov_type
    new_decorrelate_out <- structure(decorrelate_out, class = "decorrelate_list")
    return(new_decorrelate_out)
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

  # store NA values in newdata
  # all training_index and test_index values correspond to the
  # SUBSETTED data and may not correspond to the original data order
  na_index <- is.na(data[[all.vars(formula)[1]]])
  # store observed index
  if (any(na_index)) {
    newdata <- data[na_index, , drop = FALSE]
    data <- data[!na_index, , drop = FALSE]
  } else {
    newdata <- NULL
  }

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
    newdata = newdata,
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
