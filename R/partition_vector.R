#' Create a partiition vector
#'
#' @param partition_factor A partition factor (formula)
#' @param data data
#' @param newdata newdata (for prediction)
#'
#' @return
#'
#' @noRd
partition_vector <- function(partition_factor, data, newdata, reform_bar2 = NULL, partition_index_data = NULL) {
  if (is.null(partition_factor)) {
    t_partition_index <- NULL
  } else {
    if (is.null(reform_bar2)) {
      partition_factor_val <- get_randcov_name(labels(terms(partition_factor)))
      bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
      reform_bar2 <- reformulate(paste0("as.numeric(", bar_split[[2]], ")"), intercept = FALSE)
    }
    if (is.null(partition_index_data)) {
      partition_index_data <- as.vector(model.matrix(reform_bar2, data))
    }
    partition_index_newdata <- as.vector(model.matrix(reform_bar2, newdata))
    partition_index <- vapply(partition_index_newdata, function(x) ifelse(x == partition_index_data, 1, 0), numeric(length(partition_index_data)))
    t_partition_index <- Matrix(t(partition_index), sparse = TRUE)
  }
  t_partition_index
}
