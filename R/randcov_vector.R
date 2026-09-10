# prediction-time counterpart to get_randcov_list.R's fitting-time random
# effect design-matrix construction; naming conventions between the two were
# reconciled in a recent bug fix (e.g. group labels use the same
# "varname + level" convention on both sides)
#' Build a random-effect slope value vector for \code{newdata}, erroring on NA
#'
#' \code{model.matrix()} on a plain data frame (rather than a pre-built model
#' frame) defaults to \code{na.action = na.omit}, silently dropping NA rows
#' instead of keeping them -- for a single-row \code{newdata} that means a
#' 0-length result, which crashes the \code{se.fit} computation
#' (\code{randcov_newvar()}) deep inside \code{vapply()} with an uninformative
#' "result is length 0" error, and otherwise (no \code{se.fit}) silently
#' produces \code{NA} via out-of-bounds indexing (\code{get_randcov_vectors()})
#' with no warning at all. Pre-building the model frame with
#' \code{na.action = na.pass} keeps the row (as an \code{NA} value) instead of
#' dropping it, so the explicit check below can catch it with a message
#' matching the existing fixed-effect NA check in
#' \code{get_prediction_object_splm()}/\code{_spglm()}/\code{predict_block_splm()}.
#'
#' @param reform_bar1 The random effect's slope formula (\code{~ x - 1})
#' @param newdata The newdata to build the slope value from
#'
#' @return A numeric vector, one value per row of \code{newdata}
#'
#' @noRd
get_randcov_slope_val_newdata <- function(reform_bar1, newdata) {
  slope_val_newdata <- as.vector(model.matrix(reform_bar1, model.frame(reform_bar1, newdata, na.action = na.pass)))
  if (anyNA(slope_val_newdata)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }
  slope_val_newdata
}

#' Create a random effects covariance vector
#'
#' @param randcov_params A \code{cov_params} object
#' @param data Data
#' @param newdata Newdata (used for prediction)
#' @param randcov_terms An optional list (named by random effect term) of
#'   precomputed per-term context built by \code{get_extra_randcov_list()}:
#'   \code{reform_bar2}, \code{reform_bar1}, \code{group_label}, \code{xlev},
#'   \code{level_index_map}, and \code{slope_val}. When omitted, everything is
#'   derived from \code{data} directly.
#'
#' @return A random effects covariance vector
#'
#' @noRd
randcov_vector <- function(randcov_params = NULL, data, newdata, randcov_terms = NULL) {
  if (is.null(randcov_params)) {
    randcov_vectors <- NULL
  } else {
    randcov_names <- names(randcov_params)
    randcov_vectors <- lapply(
      randcov_names, get_randcov_vectors, randcov_params, data, newdata, randcov_terms
    )
    randcov_vectors <- Reduce("+", randcov_vectors)
  }
  randcov_vectors
}

#' Create the random effects covariance vector for a single random effect term
#'
#' @param randcov_name The name of a single random effect term
#' @param randcov_params A \code{randcov_params} object
#' @param data Data
#' @param newdata Newdata (used for prediction)
#' @param randcov_terms An optional list (named by random effect term) of
#'   precomputed per-term context built by \code{get_extra_randcov_list()}
#'
#' @return A sparse matrix (\code{newdata} rows by \code{data} rows) of this
#'   term's contribution to the covariance between new and observed data
#'
#' @noRd
get_randcov_vectors <- function(randcov_name, randcov_params, data, newdata, randcov_terms) {
  randcov_param <- randcov_params[randcov_name]
  bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
  term <- randcov_terms[[randcov_name]]

  if (is.null(term)) {
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
    Z_index_data_mf <- model.frame(reform_bar2, data)
    group_label <- model_matrix_group_labels(reform_bar2, data)
    xlev <- .getXlevels(terms(Z_index_data_mf), Z_index_data_mf)
    xlev_full <- .getXlevels(terms(Z_index_data_mf), rbind(Z_index_data_mf, model.frame(reform_bar2, newdata)))
    if (!identical(xlev, xlev_full)) {
      xlev <- xlev_full
    }
    # computed fresh here (once per call) when no cached context is passed in
    level_index_map <- split(seq_along(group_label), group_label)
    if (bar_split[[1]] != "1") {
      reform_bar1 <- reformulate(bar_split[[1]], intercept = FALSE)
      slope_val <- as.vector(model.matrix(reform_bar1, data))
    } else {
      reform_bar1 <- NULL
      slope_val <- NULL
    }
  } else {
    reform_bar2 <- term$reform_bar2
    group_label <- term$group_label
    xlev <- term$xlev
    # reuse the map precomputed once in get_extra_randcov_list() instead of
    # rescanning all of data on every call (e.g., once per prediction row)
    level_index_map <- term$level_index_map
    if (is.null(level_index_map)) {
      level_index_map <- split(seq_along(group_label), group_label)
    }
    reform_bar1 <- term$reform_bar1
    slope_val <- term$slope_val
  }

  group_label_newdata <- model_matrix_group_labels(reform_bar2, newdata, xlev = xlev, na_pass = TRUE)

  # for each newdata row, look up the (typically few) data rows sharing its group
  # label via a hash/list lookup instead of comparing against every data row --
  # this is the change that scales with match count rather than nrow(data)
  n_obs <- length(group_label)
  n_new <- length(group_label_newdata)
  matches <- vector("list", n_new)
  non_na_new <- !is.na(group_label_newdata)
  matches[non_na_new] <- level_index_map[group_label_newdata[non_na_new]]
  match_lengths <- lengths(matches)
  obs_idx <- unlist(matches, use.names = FALSE)
  new_idx <- rep(seq_len(n_new), match_lengths)

  if (!is.null(reform_bar1)) {
    slope_val_newdata <- get_randcov_slope_val_newdata(reform_bar1, newdata)
    # equivalent to the old dense Z_index (cov zeroed off-group) swept by
    # slope_val_newdata and slope_val, but only ever computed at the nonzero
    # (matched) positions
    x_vals <- randcov_param * slope_val[obs_idx] * slope_val_newdata[new_idx]
  } else {
    x_vals <- rep(randcov_param, length(obs_idx))
  }

  Matrix::sparseMatrix(i = new_idx, j = obs_idx, x = x_vals, dims = c(n_new, n_obs))
}

#' The marginal variance random effects contribute to a single new observation
#'
#' @param randcov_params A \code{cov_params} object
#' @param newdata A single-row data frame for the new observation being predicted
#' @param randcov_terms An optional list (named by random effect term) of
#'   precomputed per-term context built by \code{get_extra_randcov_list()}
#'
#' @details For a random intercept, \code{Var(b_g)} is just the term's variance
#'   parameter. For a random slope (\code{x | group}), \code{Var(b_g * x0) =
#'   sigma^2 * x0^2} depends on the new observation's covariate value, so it
#'   cannot be read off \code{randcov_params} directly the way it can for an
#'   intercept -- unlike the training-data diagonal (which factors through
#'   \code{ZZt} and so already reflects each observation's covariate value),
#'   nothing upstream of this computes that value for a location that isn't in
#'   the data.
#'
#' @return A single numeric value: the total random effect contribution to
#'   \code{Var(Y0)} for this new observation
#'
#' @noRd
randcov_newvar <- function(randcov_params = NULL, newdata, randcov_terms = NULL) {
  if (is.null(randcov_params)) {
    return(0)
  }
  randcov_names <- names(randcov_params)
  vars <- vapply(randcov_names, function(randcov_name) {
    randcov_param <- as.numeric(randcov_params[randcov_name])
    bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
    if (bar_split[[1]] == "1") {
      return(randcov_param)
    }
    term <- randcov_terms[[randcov_name]]
    reform_bar1 <- if (is.null(term)) reformulate(bar_split[[1]], intercept = FALSE) else term$reform_bar1
    slope_val_newdata <- get_randcov_slope_val_newdata(reform_bar1, newdata)
    randcov_param * slope_val_newdata^2
  }, numeric(1))
  sum(vars)
}
