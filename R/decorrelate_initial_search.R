#' Evaluate a decorrelation parameter grid on one training/test split
#'
#' The core of \code{\link{decorrelate}()}'s grid search: for a single
#' training/test split (one element of the resolved \code{training} list;
#' see \code{\link{get_training_list}()}), decorrelates the training data
#' under every candidate parameter set in \code{grid}, fits the machine
#' learning algorithm to each, decorrelates+predicts the test data, and
#' recorrelates the predictions to compute test-set fit statistics (bias,
#' MSPE, RMSPE, cor2) for that parameter set. \code{\link{decorrelate}()}
#' calls this once per training replicate/fold and averages the results
#' across calls when there is more than one.
#'
#' @param formula,data,spcov_type,spcov_params,algorithm,statistic,anisotropy,random,randcov_params,partition_factor,ordering,local,grid,dense_grid
#'   See \code{\link{decorrelate}()}.
#' @param training_list A single \code{list(training_index, test_index)}
#'   element from a resolved \code{\link{get_training_list}()} result.
#' @param add_iid Whether to append an untransformed baseline row to the
#'   grid when one is not already supplied; see \code{\link{decorrelate}()}.
#' @param ... Additional arguments passed to the machine learning algorithm.
#'
#' @return A list with elements \code{params_list} (the grid rows as
#'   \code{spcov_params}/\code{randcov_params} objects), \code{grid} (the
#'   grid with \code{bias}/\code{MSPE}/\code{RMSPE}/\code{cor2} columns
#'   appended), and \code{training} (\code{training_list}, passed through).
#'
#' @noRd
decorrelate_initial_search <- function(formula, data, spcov_type, spcov_params, xcoord, ycoord, algorithm, statistic, training_list, anisotropy, random, randcov_params, partition_factor, ordering = "maxmin", local, grid, dense_grid, add_iid, ...) {


  data_training <- data[training_list$training_index, , drop = FALSE]
  data_test <- data[training_list$test_index, , drop = FALSE]
  yval <- model.response(model.frame(formula, data = data_test))

  # build a default candidate grid unless the user supplied one; a supplied
  # grid is validated and column-matched against what decorrelate_grid_internal()
  # would have produced, so downstream code can treat both cases identically
  grid_compare <- decorrelate_grid_internal(
    formula = formula,
    data = data,
    spcov_type = spcov_type,
    spcov_params = spcov_params,
    xcoord = xcoord,
    ycoord = ycoord,
    anisotropy = anisotropy,
    random = random,
    randcov_params = randcov_params,
    dense_grid = dense_grid,
    add_iid = add_iid,
    warn = FALSE
  )
  if (is.null(grid)) {
    grid <- grid_compare
  } else {
    check_grid_legal(grid, random)
    if (!"rotate" %in% names(grid)) grid$rotate <- 0
    if (!"scale" %in% names(grid)) grid$scale <- 1
    if (any(! names(grid_compare) %in% names(grid))) stop("Invalid grid column names. Column names in grid must contain all column names returned by decorrelate_grid(formula, data, ...).", call. = FALSE)
    grid <- grid[, names(grid) %in% names(grid_compare), drop = FALSE]
  }

  params_list <- get_params_list(grid, random, randcov_params)

  # the expensive setup shared by every grid row (data object, ordering,
  # local approximation) is computed once here via part1 and reused by
  # part2 inside the lapply() below -- see decorrelate_data_internal_part1()
  decorrelate_part1 <- decorrelate_data_internal_part1(
    formula = formula,
    data = data_training,
    spcov_type = spcov_type,
    xcoord = xcoord,
    ycoord = ycoord,
    random = random,
    partition_factor = partition_factor,
    ordering = ordering,
    local = local,
    ...
  )

  # for each candidate parameter set: decorrelate the training data, fit the
  # machine learning algorithm, decorrelate the test data the same way, predict, then
  # recorrelate the predictions back to the original (response) scale before
  # comparing to the held-out yval -- fit statistics are always computed on
  # the recorrelated (original-scale) predictions, not the decorrelated ones
  out <- lapply(params_list, function(x) {
    tdata_training <- decorrelate_data_internal_part2(
      spcov_params = x$spcov_params,
      randcov_params = x$randcov_params,
      decorrelate_part1_object = decorrelate_part1,
      ...
    )

    fit <- fit_decorrelate_algorithm(tdata_training, algorithm, ...)
    tdata_test <- decorrelate_newdata(tdata_training, newdata = data_test)
    preds <- predict_decorrelate_algorithm(fit, tdata_test, algorithm, ...)
    sp_decorr_preds <- recorrelate_newdata(tdata_test, preds)
    errors <- yval - sp_decorr_preds
    bias <- mean(errors)
    MSPE <- mean(errors^2)
    RMSPE <- sqrt(MSPE)
    cor2 <- suppressWarnings(cor(yval, sp_decorr_preds))^2 # warning for iid data when sp_decorr_preds = 0
    list(bias = bias, MSPE = MSPE, RMSPE = RMSPE, cor2 = cor2)
  })
  grid$bias <- unlist(lapply(out, function(x) x$bias))
  grid$MSPE <- unlist(lapply(out, function(x) x$MSPE))
  grid$RMSPE <- unlist(lapply(out, function(x) x$RMSPE))
  grid$cor2 <- unlist(lapply(out, function(x) x$cor2))

  list(params_list = params_list,
       grid = grid, training = training_list
      )
}

#' Convert a decorrelation parameter grid's rows into parameter objects
#'
#' Each row of a \code{\link{decorrelate_grid}()}-style grid stores spatial
#' covariance and random effect variance values as plain numeric columns;
#' this converts every row into the \code{\link{spcov_params}()} object (and,
#' if random effects are present, named variance vector) that
#' \code{\link{decorrelate_data_internal_part2}()} and friends expect.
#'
#' @param grid A grid as returned by \code{\link{decorrelate_grid}()}/
#'   \code{\link{decorrelate_grid_internal}()}.
#' @param random The \code{random} formula, or \code{NULL}.
#' @param randcov_params A \code{randcov_params} object/vector, or \code{NULL}
#'   (only used to detect whether random effects are in play; per-row values
#'   come from \code{grid} itself).
#'
#' @return A list (one element per grid row) of
#'   \code{list(spcov_params, randcov_params)}.
#'
#' @noRd
get_params_list <- function(grid, random, randcov_params) {
  params_list <- lapply(seq(1, NROW(grid)), function(x) {
    x <- grid[x, ]
    spcov_type <- x[["spcov_type"]]
    if ("extra" %in% names(x) && spcov_type != "none") {
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
      remove_cols <- c("spcov_type", "de", "ie", "range", "rotate", "scale", "bias", "MSPE", "RMSPE", "cor2")
      if ("extra" %in% names(x)) remove_cols <- c(remove_cols, "extra")
      randcov_params_val <- unlist(x[, -which(names(x) %in% remove_cols), drop = FALSE])
      names(randcov_params_val) <- paste("", names(randcov_params_val), "", sep = "")
    } else {
      randcov_params_val <- NULL
    }
    list(spcov_params = spcov_params_val, randcov_params = randcov_params_val)
  })
}
