#' Apply the Spatial Decorrelation Transformation for Machine Learning Models
#'
#' @description Apply the spatial decorrelation transformation for
#'   point-referenced data, allowing for random effects,
#'   anisotropy, partition factors, and big data methods.
#'
#' @param formula A two-sided linear formula describing the fixed effect structure
#'   of the model, with the response to the left of the \code{~} operator and
#'   the terms on the right, separated by \code{+} operators. \code{.} on the
#'   right-hand side represents every variable in \code{data} except the
#'   response and the x-coordinate/y-coordinate columns (\code{xcoord}/\code{ycoord},
#'   or, for an \code{sf} object, the geometry column), which are never
#'   included via \code{.} (though they may still be given explicitly).
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
#' @param statistic The statistic used to evaluate fit in the test data. Available options
#'   include \code{"bias"} (mean bias), \code{"MSPE"} (mean-squared-prediction error),
#'    \code{"RMSPE"} (root-mean-squared-prediction error)
#'   and \code{"cor2"} (the predictive R-squared; i.e., the
#'   squared correlation between observations and predictions).
#' @param training An list controlling how the training and test data are assigned
#'   when evaluating test data performance.
#'   The following arguments detail this process:
#'   \itemize{
#'    \item \code{method}: The method used to evaluate test data performance.
#'      \code{"split"} will split \code{data} up into distinct training and test sets
#'      proportionally based on \code{p}. \code{"cv"} will split \code{data} up
#'      via k-fold cross validation based on \code{folds}, the number of folds.
#'    \item \code{p}: The proportion (a numeric vector between zero and one) of observations in \code{data} that should
#'      be assigned to the training data. The default is 0.8, which means that
#'      80% of the observations are assigned to the training data and 20% to the
#'      test data. Ignored if \code{training_index} or \code{test_index} are provided.
#'    \item \code{replicate}: The number of times to replicate \code{"split"} with different random training and test assignments.
#'    \item \code{folds}: The number of folds to use in cross-validation. Requires \code{method = "cv"}. Ignored if \code{folds_index} is specified. The default is 5, matching the 80/20 default split above.
#'    \item \code{training_index}: A numeric vector that specifies which rows (i.e., indices)
#'      of \code{data} should be assigned to the training data. If omitted, defaults
#'      to the rows which are not already included in \code{test_index}.
#'    \item \code{test_index}: A numeric vector that specifies which rows (i.e., indices)
#'      of \code{data} should be assigned to the test data. If omitted, defaults
#'      to the rows which are not already included in \code{training_index}.
#'    \item \code{folds_index}: A numeric vector that specifies which rows
#'      of \code{data} are associated with each cross-validation fold. Requires
#'      \code{method = "cv"}.
#'   }
#'   If omitted, \code{training} is transformed into
#'   \code{list(method = "split", p = 0.8, replicate = 1)}.
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
#'   include \code{"grts"}, \code{"maxmin"}, \code{"middleout"},
#'   \code{"outsidein"}, \code{"coordinate"}, \code{"random"}, and \code{"none"}.
#'   \code{"grts"} applies ordering using a spatially balanced GRTS sample via \code{spsurvey::grts()}.
#'   \code{"maxmin"} applies maximum minimum distance ordering via \code{GPvecchia::order_maxmin_exact()}.
#'   \code{"middleout"} applies middle out ordering via \code{GPvecchia::order_middleout()}.
#'   \code{"outsidein"} applies middle out ordering via \code{GPvecchia::order_outsidein()}.
#'   \code{"coordinate"} applies middle out ordering via \code{GPvecchia::order_coordinate(..., coordinate = c(1, 2))},
#'   which orders from bottom-left to top-right of the spatial domain.
#'   \code{"random"} applies a completely random ordering.
#'   \code{"none"} applies no random ordering.
#'   The default is \code{"maxmin"} unless there are multiple observations at a single
#'   location, in which case the default is \code{"grts"}.
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
#' @param grid An explicit grid of parameter values by which to evaluate fit. The
#'   names of \code{grid} must contain all the names returned by \code{decorrelate_grid(formula, data, ...)}.
#' @param dense_grid If \code{grid} is not provided, \code{dense_grid} is a logical
#'   which controls the density of the constructed grid to be evaluated. If
#'   \code{dense_grid} is \code{TRUE}, a denser grid is used. If \code{dense_grid}
#'   is \code{FALSE}, a sparser grid is used. By default, \code{dense_grid}
#'   is \code{FALSE} when the sample size is greater than 5,000 and \code{TRUE}
#'   otherwise.
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
#'     \item matern: \eqn{2^{1 - extra}/ \Gamma(extra) * \alpha^{extra} * Bk(\alpha, extra)}, \eqn{\alpha = (2extra)^{0.5} * \eta}, Bk is Bessel-K function with order \eqn{1/5 \le extra \le 5}
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
#'   \code{training} Details: When \code{replicate} or the number of cross validation folds is at
#'     least two, there are separate grids evaluated for each replication (or fold). Statistics in each grid are
#'     averaged across replications (or folds) to determine a final grid ranked by \code{statistic}.
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
#'     \item \code{algorithm}: The machine learning algorithm used.
#'     \item \code{decorrelate_data}: The output of [decorrelate_data()] applied to \code{data}.
#'     \item \code{fit}: The fitted machine learning model object applied to the decorrelated data.
#'     \item \code{grid}: If used, the grid of spatial decorrelation parameters evaluated and their corresponding
#'       metrics when applied to the test data.
#'     \item \code{newdata}: The rows of \code{data} that have \code{NA} response values and are stored as prediction data.
#'     \item \code{training}: If used, the observations assigned to each training and test data set.
#'     \item \code{test}: If used, a list with the lowest (absolute) mean bias (bias), mean-squared-prediction error (MSPE),
#'       root-mean-squared-prediction error (RMSPE), and predictive R-squared (cor2).
#'   }
#'
#' @name decorrelate
#' @order 1
#' @export
#'
#' @references Matthew J. Heaton, Andrew Millane, and Jake S. Rhodes. 2025. A Scalable
#'   Spatial Decorrelation Preprocessing Approach for Machine and Deep Learning.
#'   \emph{Journal of Data Science}. 1-15, DOI 10.6339/25-JDS1210
#'
#' @examples
#' decorr <- decorrelate(log_cond ~ temp, data = lake, spcov_type = "exponential")
#' tidy(decorr$grid)
decorrelate <- function(formula, data, spcov_type, spcov_params, xcoord, ycoord, algorithm = "ranger", statistic = "RMSPE", training, evaluate_test, anisotropy = FALSE, random, randcov_params, partition_factor, ordering, local, grid, dense_grid, ...) {

  # check default grid
  if (!missing(spcov_type) && !missing(grid) && missing(spcov_params)) {
    message("Both spcov_type and grid provided. grid overriding spcov_type.")
  }

  if (missing(random)) random <- NULL
  if (missing(randcov_params)) randcov_params <- NULL
  if (!is.null(randcov_params) && is.null(random)) {
    random <- reformulate(names(randcov_params))
  }

  # set exponential as default if nothing specified
  if (missing(spcov_type) && missing(spcov_params)) {
    if (missing(grid)) {
      spcov_type <- "exponential"
      message("No spatial covariance type provided. Assuming \"exponential\".")
    } else {
      # infer spcov_type from a user-supplied grid instead: a grid with only
      # "none" is a single type; a grid mixing "none" with one other type
      # represents that other type plus an explicit untransformed baseline
      # row (de/range are filled in as the "no transformation" placeholder
      # values so the grid stays legal for the estimation code below)
      check_grid_legal(grid, random)
      unq_spcov_type <- unique(grid$spcov_type)
      if (length(unq_spcov_type) == 1) {
        spcov_type <- unq_spcov_type
      } else if (length(unq_spcov_type) == 2) {
        grid$de[grid$spcov_type == "none"] <- 0
        grid$range[grid$spcov_type == "none"] <- Inf
        spcov_type <- unq_spcov_type[unq_spcov_type != "none"]
      } else {
        stop("No rows in grid.", call. = FALSE)
      }
    }
  }

  if (! statistic %in% c("bias", "MSPE", "RMSPE", "cor2")) stop("statistic must be \"bias\", \"MSPE\", \"RMSPE\" or \"cor2\".", call. = FALSE)

  # multiple spcov_type values: recurse once per type (reusing the same
  # resolved training split across all of them, so results are comparable),
  # collecting the fits into a named "decorrelate_list" instead of running
  # the rest of this function once per type
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

  xcoord <- if (missing(xcoord)) NULL else as.character(substitute(xcoord))
  ycoord <- if (missing(ycoord)) NULL else as.character(substitute(ycoord))
  if (missing(spcov_type) && !missing(spcov_params)) {
    spcov_type <- class(spcov_params)
  }
  if (missing(spcov_params)) spcov_params <- NULL
  if (missing(partition_factor)) partition_factor <- NULL
  if (missing(ordering)) ordering <- NULL
  if (missing(local)) local <- NULL
  if (missing(training)) training <- NULL
  if (missing(evaluate_test)) evaluate_test <- FALSE
  if (missing(grid)) grid <- NULL

  # write a decorrelate_checks()?

  # store NA values in newdata
  # all training_index and test_index values correspond to the
  # SUBSETTED data and may not correspond to the original data order
  # rows with a missing response are treated as prediction data (mirroring
  # splm()/spglm()'s newdata handling), not as part of the estimation sample
  na_index <- is.na(data[[all.vars(formula)[1]]])
  # store observed index
  if (any(na_index)) {
    newdata <- data[na_index, , drop = FALSE]
    data <- data[!na_index, , drop = FALSE]
  } else {
    newdata <- NULL
  }
  training <- get_training_list(training, data)

  if (missing(dense_grid)) {
    if (NROW(data) <= 5000) {
      dense_grid <- TRUE
    } else {
      dense_grid <- FALSE
    }
  }

  # a grid search is only needed when at least one decorrelation parameter
  # (spatial and/or random effect) is not already fully specified by the user
  if (is.null(spcov_params) || (!is.null(random) && is.null(randcov_params))) {
    evaluate_test <- TRUE
  }
  if (evaluate_test) {

    # add_iid: whether to append a fully "untransformed" (spcov_type =
    # "none", ie = 1) baseline row to the grid so the transform can be
    # compared against applying the machine learning algorithm to the raw data. Only added
    # when parameters are actually being estimated (an explicit grid or
    # fully-known spcov_params/randcov_params means the user already chose
    # what to evaluate, so this default baseline is skipped)
    if (is.null(random)) {
      add_iid <- ifelse(!is.null(spcov_params) || !is.null(grid), FALSE, TRUE)
    } else {
      add_iid <- ifelse((!is.null(spcov_params) && !is.null(randcov_params)) || !is.null(grid), FALSE, TRUE)
    }

    # evaluate the candidate grid once per training/test split (replicated
    # "split" or "cv" folds), then combine results below
    names_training <- names(training$training)
    init <- lapply(names_training, function(x) {
      decorrelate_initial_search(
        formula = formula,
        data = data,
        spcov_type = spcov_type,
        spcov_params = spcov_params,
        xcoord,
        ycoord,
        algorithm = algorithm,
        statistic = statistic,
        training_list = training$training[[x]],
        anisotropy = anisotropy,
        random = random,
        randcov_params = randcov_params,
        partition_factor = partition_factor,
        ordering = ordering,
        local = local,
        grid = grid,
        dense_grid = dense_grid,
        add_iid = add_iid,
        ...
      )
    })
    if (length(names_training) == 1) {
      grid <- init[[1]]$grid
    } else {
      # multiple splits/folds: average each parameter set's fit statistics
      # across replicates before picking a winner, rather than picking a
      # winner per replicate and averaging afterward
      grids <- do.call(rbind, lapply(init, function(x) x$grid))
      form <- cbind(bias, MSPE, RMSPE, cor2) ~ .
      grid <- aggregate(x = form, data = grids, FUN = mean)
      # grid[, "RMSPE"] <- sqrt(grid[, "MSPE"])
    }

    # pick the best-fitting parameter set: highest cor2, or smallest
    # (absolute, for bias) error otherwise
    if (statistic == "cor2") {
      best_val <- which.max(grid[[statistic]])
    } else if (statistic == "bias") {
      best_val <- which.min(abs(grid[[statistic]]))
    } else {
      best_val <- which.min(grid[[statistic]])
    }

    test <- list(
      bias = grid$bias[best_val],
      MSPE = grid$MSPE[best_val],
      RMSPE = grid$RMSPE[best_val],
      cor2 = grid$cor2[best_val],
      statistic = statistic
    )

    # convert the winning grid row back into spcov_params/randcov_params
    # objects usable by decorrelate_data_internal() below
    params_list <- get_params_list(grid, random, randcov_params)
    spcov_params <- params_list[[best_val]]$spcov_params
    randcov_params <- params_list[[best_val]]$randcov_params

    if (statistic == "cor2") {
      grid <- grid[order(grid[[statistic]], decreasing = TRUE), , drop = FALSE]
    } else if (statistic == "bias") {
      grid <- grid[order(abs(grid[[statistic]])), , drop = FALSE]
    } else {
      grid <- grid[order(grid[[statistic]]), , drop = FALSE]
    }
    row.names(grid) <- NULL
  } else {
    # spcov_params (and randcov_params, if relevant) were fully specified by
    # the user and evaluate_test was left FALSE, so there is nothing to
    # search or report -- just move on to fitting with the known parameters
    grid <- NULL
    training <- NULL
    test <- NULL
  }

  decorr <- decorrelate_data_internal(
    formula = formula,
    data = data,
    spcov_params = spcov_params,
    xcoord = xcoord,
    ycoord = ycoord,
    random = random,
    randcov_params = randcov_params,
    partition_factor = partition_factor,
    ordering = ordering,
    local = local,
    ...
  )

  fit <- fit_decorrelate_algorithm(decorr, algorithm, ...)

  # set grid class for printing later
  if (!is.null(grid)) {
    grid <- structure(grid, class = c("decorrelate_grid", class(grid)), statistic = statistic)
  }

  obj <- list(
    algorithm = algorithm,
    call = match.call(),
    decorrelate_data = decorr,
    fit = fit,
    grid = grid,
    newdata = newdata,
    training = training,
    test = test
  )
  new_obj <- structure(obj, class = "decorrelate")
  new_obj
}

#' @rdname decorrelate
#' @param x An object from \code{object$grid}.
#' @param sort_by Sort by a specific row in \code{x}. Fit statistics are
#'   \code{"bias"}, \code{"MSPE"}, \code{"RMSPE"}, and \code{"cor2"}.
#'   The default is \code{"MSPE"}.
#' @param decreasing Whether \code{sort_by} should sort by decreasing order? If
#' \code{sort_by = "cor2"}, the default is \code{TRUE}; otherwise it is \code{FALSE}.
#' @method tidy decorrelate_grid
#' @order 2
#' @export
tidy.decorrelate_grid <- function(x, sort_by, decreasing, ...) {

  if (missing(sort_by)) {
    sort_by <- attr(x, "statistic")
  }

  if (!sort_by %in% names(x)) {
    stop("sort_by must be a variable in x.", call. = FALSE)
  }

  if (missing(decreasing)) {
    decreasing <- ifelse(sort_by == "cor2", TRUE, FALSE)
  } else {
    if (!is.logical(decreasing)) {
      stop("decreasing must be TRUE or FALSE.", call. = FALSE)
    }
  }


  x <- x[order(x[[sort_by]], decreasing = decreasing), , drop = FALSE]

  # identify the "no transformation" baseline row(s) added by
  # decorrelate()'s add_iid logic: spcov_type = "none" with a trivial nugget
  # (ie = 1), infinite range, no anisotropy, and (if present) every random
  # effect variance at zero -- i.e. every parameter that would otherwise
  # decorrelate the data is switched off
  untransformed_index <- x$spcov_type == "none" & x$de == 0 & x$ie == 1 & x$range == Inf & x$rotate == 0 & x$scale == 1
  if ("extra" %in% names(x)) {
    untransformed_index <- untransformed_index & x$extra == 0
  }
  standard_names <- c("spcov_type", "de", "ie", "range", "extra", "rotate", "scale", "bias", "MSPE", "RMSPE", "cor2")
  random_names <- names(x)[!names(x) %in% standard_names]
  if (length(random_names) > 0) {
    random_index <- apply(do.call(cbind, lapply(random_names, function(y) x[[y]] == 0)), 1, all)
    untransformed_index <- untransformed_index & random_index
  }

  if (all(x$rotate == 0 & x$scale == 1)) {
    x$rotate <- NULL
    x$scale <- NULL
  }

  # relabel those rows for display: literal parameter values like ie = 1 are
  # misleading for a row that applies no transformation at all, so report
  # NA for every decorrelation parameter instead
  if (any(untransformed_index)) {
    x[untransformed_index, "spcov_type"] <- "no transformation"
    x[untransformed_index, "de"] <- NA
    x[untransformed_index, "ie"] <- NA
    x[untransformed_index, "range"] <- NA
    if ("extra" %in% names(x)) {
      x[untransformed_index, "extra"] <- NA
    }
    if (length(random_names) > 0) {
      for (y in random_names) {
        x[untransformed_index, y] <- NA
      }
    }
    if (any(c("rotate", "scale") %in% names(x))) {
      x[untransformed_index, "rotate"] <- NA
      x[untransformed_index, "scale"] <- NA
    }
  }
  tibble::tibble(x)
}


#' Fit a machine learning algorithm to spatially decorrelated data
#'
#' Thin dispatcher that fits \code{decorrelate_data$tX}/\code{ty} (the
#' spatially decorrelated design matrix and response; see
#' \code{\link{decorrelate_data}()}) using whichever machine learning package
#' \code{algorithm} names, after checking that package is installed.
#'
#' @param decorrelate_data A \code{\link{decorrelate_data}()} object.
#' @param algorithm \code{"ranger"}, \code{"randomForest"}, or \code{"xgboost"}.
#' @param ... Additional arguments passed to the underlying fitting function.
#'
#' @return The fitted model object returned by the chosen algorithm's package.
#'
#' @noRd
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

#' Predict from a fitted machine learning algorithm on decorrelated test data
#'
#' Counterpart to \code{\link{fit_decorrelate_algorithm}()} used while
#' evaluating the parameter grid: each package's \code{predict()} method
#' takes its transformed-newdata argument under a different name
#' (\code{data} for ranger, \code{newdata} for randomForest/xgboost), so this
#' normalizes the call. Predictions here are still on the decorrelated scale;
#' the caller recorrelates them before computing test-set fit statistics.
#'
#' @param decorrelate_data A fitted model object from
#'   \code{\link{fit_decorrelate_algorithm}()}.
#' @param tdata_test A list with element \code{tX_newdata}, the spatially
#'   decorrelated design matrix for the test data.
#' @param algorithm \code{"ranger"}, \code{"randomForest"}, or \code{"xgboost"}.
#' @param ... Additional arguments passed to the underlying \code{predict()}.
#'
#' @return A numeric vector of decorrelated-scale predictions for the test data.
#'
#' @noRd
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

  # if (algorithm == "nnet") {
  #   preds <- predict(decorrelate_data, newdata = tdata_test$tX_newdata, ...)
  # }

  preds
}

#' Build the resolved training/test split(s) used to evaluate the
#' decorrelation parameter grid
#'
#' Fills in defaults for (and validates) the \code{training} argument to
#' \code{\link{decorrelate}()}, normalizing every supported specification
#' (a train proportion, explicit \code{training_index}/\code{test_index}, a
#' replicated random split, or k-fold cross-validation) down to a common
#' \code{training$training} list of \code{list(training_index, test_index)}
#' elements -- one element per replicate/fold -- that
#' \code{\link{decorrelate_initial_search}()} iterates over.
#'
#' @param training A list; see the \code{training} argument to
#'   \code{\link{decorrelate}()}.
#' @param data The (non-missing-response) data being split.
#'
#' @return The resolved \code{training} list, with a \code{training$training}
#'   element as described above.
#'
#' @noRd
get_training_list <- function(training, data) {

  if (is.null(training)) {
    training <- list()
  }

  if (!is.list(training)) {
    stop("training must be a list")
  }

  names_training <- names(training)

  if (!"method" %in% names_training) {
    training$method <- "split"
  }

  if (training$method == "split") {
    # do stuff

    if (!"p" %in% names_training) {
      training$p <- 0.8
    }

    if (!"replicate" %in% names_training) {
      training$replicate <- 1
    }

    # three ways to specify a split: a train proportion p (drawing
    # `replicate` independent random splits), an explicit training_index
    # (test = everything else, one split only), or an explicit test_index
    # (training = everything else, one split only)
    if (!"training_index" %in% names_training && !"test_index" %in% names_training) {
      n <- NROW(data)
      index <- seq(1, n)
      n_training <- floor(n * training$p)
      n_test <- n - n_training
      training$training <- lapply(seq(1, training$replicate), function(x) {
        index <- sample(index)
        training_index <- index[seq(1, n_training)]
        test_index <- index[seq(n_training + 1, n)]
        list(training_index = training_index, test_index = test_index)
      })
    } else if ("training_index" %in% names_training && !"test_index" %in% names_training) {
      n <- NROW(data)
      index <- seq(1, n)
      training_index <- training$training_index
      test_index <- index[-training$training_index]
      training$training <- list(list(training_index = training_index, test_index = test_index))
      training$training_index <- NULL
    } else if (!"training_index" %in% names_training && "test_index" %in% names_training) {
      n <- NROW(data)
      index <- seq(1, n)
      test_index <- training$test_index
      training_index <- index[-training$test_index]
      training$training <- list(list(training_index = training_index, test_index = test_index))
      training$test_index <- NULL
    }

    # sort
    training$training <- lapply(training$training, function(x) {
      training_index <- sort(x$training_index)
      test_index <- sort(x$test_index)
      if (length(training_index) == 0) stop("No observations detected in training data.", call. = FALSE)
      if (length(test_index) == 0) stop("No observations detected in test data.", call. = FALSE)
      if (any(duplicated(training_index))) stop("Cannot have duplicated rows in training_index.", call. = FALSE)
      if (any(duplicated(test_index))) stop("Cannot have duplicated rows in test_index.", call. = FALSE)
      list(training_index = training_index, test_index = test_index)
    })
    names(training$training) <- as.character(seq(1, training$replicate))

  }

  if (training$method == "cv") {
    # do stuff

    n <- NROW(data)
    index <- seq(1, n)

    if (!"folds" %in% names_training) {
      training$folds <- 5
    }

    if (!"folds_index" %in% names_training) {
      folds_index <- rep(seq(1, training$folds), length.out = n)
      training$folds_index <- sample(folds_index)
    }

    # each fold takes a turn as the test set, with the remaining folds as
    # training -- same list-of-splits shape as the "split" method above, so
    # decorrelate_initial_search() doesn't need to know which method produced it
    unq_folds <- sort(unique(training$folds_index))
    training$training <- lapply(unq_folds, function(x) {
      in_fold <- training$folds_index == x
      training_index <- index[!in_fold]
      test_index <- index[in_fold]
      if (length(training_index) == 0) stop("No observations detected in training data.", call. = FALSE)
      if (length(test_index) == 0) stop("No observations detected in test data.", call. = FALSE)
      if (any(duplicated(training_index))) stop("Cannot have duplicated rows in training_index.", call. = FALSE)
      if (any(duplicated(test_index))) stop("Cannot have duplicated rows in test_index.", call. = FALSE)
      list(training_index = training_index, test_index = test_index)
    })
    names(training$training) <- unq_folds

  }

  training

}

#' Validate a user-supplied \code{decorrelate()} grid
#'
#' A grid may only mix spatial covariance types in one specific way: a single
#' "real" \code{spcov_type} plus an untransformed \code{"none"}/\code{ie = 1}
#' baseline row for comparison (see the \code{add_iid} logic in
#' \code{\link{decorrelate}()}). Anything else -- more than two distinct
#' types, a second type that isn't \code{"none"}, or a \code{"none"} row that
#' isn't actually the trivial baseline (\code{ie != 1}, and only relevant
#' without random effects) -- is rejected.
#'
#' @param grid A candidate grid, as returned by \code{\link{decorrelate_grid}()}
#'   or supplied directly to \code{\link{decorrelate}()}.
#' @param random The \code{random} formula, or \code{NULL}.
#'
#' @return Invisibly \code{NULL}; called for its error-checking side effect.
#'
#' @noRd
check_grid_legal <- function(grid, random) {
  unq_spcov_type <- unique(grid$spcov_type)
  if (length(unq_spcov_type) > 2 || (length(unq_spcov_type) == 2) && (!"none" %in% unq_spcov_type) || (length(unq_spcov_type) == 2) && ("none" %in% unq_spcov_type) && any(grid$ie[grid$spcov_type == "none"] != 1) && (is.null(random))) {
    stop("All spatial covariance types must be the same. The single exception is \"none\" when ie = 1 without random effects, which indicates no transformation.", call. = FALSE)
  }
}
