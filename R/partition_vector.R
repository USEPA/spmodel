#' Create a partiition vector
#'
#' @param partition_factor A partition factor (formula)
#' @param data data
#' @param newdata newdata (for prediction)
#' @param reform_bar2 An optional precomputed one-sided formula for the partition
#'   grouping variable(s)
#' @param partition_index_data An optional precomputed context built by
#'   \code{get_extra_partition_list()}: \code{group_label}, \code{xlev}, and
#'   \code{level_index_map}. When omitted, everything is derived from
#'   \code{data} directly.
#'
#' @return A partition vector for use with prediction
#'
#' @noRd
partition_vector <- function(partition_factor, data, newdata, reform_bar2 = NULL, partition_index_data = NULL) {
  # rectangular analog of partition_matrix(): an n_new x n_obs 0/1 indicator
  # of which observed rows share each newdata row's partition group, used to
  # zero out covariance between an observation to predict and any observed
  # data outside its partition
  if (is.null(partition_factor)) {
    return(NULL)
  }

  if (is.null(reform_bar2)) {
    partition_factor_val <- get_randcov_name(labels(terms(partition_factor)))
    bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
  }

  if (is.null(partition_index_data)) {
    p_index_data_mf <- model.frame(reform_bar2, data)
    p_index_data_terms <- terms(p_index_data_mf)
    group_label <- model_matrix_group_labels(reform_bar2, data)
    xlev <- .getXlevels(p_index_data_terms, p_index_data_mf)
    xlev_full <- .getXlevels(p_index_data_terms, rbind(p_index_data_mf, model.frame(reform_bar2, newdata)))
    if (!identical(xlev, xlev_full)) {
      xlev <- xlev_full
    }
    # computed fresh here (once per call) when no cached context is passed in
    level_index_map <- split(seq_along(group_label), group_label)
  } else {
    group_label <- partition_index_data$group_label
    xlev <- partition_index_data$xlev
    # reuse the map precomputed once in get_extra_partition_list() instead of
    # rescanning all of data on every call (e.g., once per prediction row)
    level_index_map <- partition_index_data$level_index_map
    if (is.null(level_index_map)) {
      level_index_map <- split(seq_along(group_label), group_label)
    }
  }

  group_label_newdata <- model_matrix_group_labels(reform_bar2, newdata, xlev = xlev, na_pass = TRUE)

  # for each newdata row, look up the (typically few) data rows sharing its
  # partition label via a hash/list lookup instead of comparing against every
  # data row -- this is the change that scales with match count rather than nrow(data)
  n_obs <- length(group_label)
  n_new <- length(group_label_newdata)
  matches <- vector("list", n_new)
  non_na_new <- !is.na(group_label_newdata)
  matches[non_na_new] <- level_index_map[group_label_newdata[non_na_new]]
  match_lengths <- lengths(matches)
  obs_idx <- unlist(matches, use.names = FALSE)
  new_idx <- rep(seq_len(n_new), match_lengths)

  # build the n_new x n_obs indicator directly from the (row, column) pairs
  # found above rather than filling a dense matrix and checking equality
  # elementwise, which would be wasteful given most entries are zero
  Matrix::sparseMatrix(i = new_idx, j = obs_idx, x = rep(1, length(obs_idx)), dims = c(n_new, n_obs))
}
