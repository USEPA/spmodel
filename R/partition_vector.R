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
      reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
    }
    if (is.null(partition_index_data)) {
      p_index_data_mf <- model.frame(reform_bar2, data)
      p_index_data_mx <- model.matrix(reform_bar2, p_index_data_mf)
      p_index_data_names <- colnames(p_index_data_mx)
      p_index_data_split <- split(p_index_data_mx, seq_len(NROW(p_index_data_mx)))
      p_index_data_xlev <- .getXlevels(terms(p_index_data_mf), p_index_data_mf)
      p_index_data_vals <- p_index_data_names[vapply(p_index_data_split, function(y) which(as.logical(y)), numeric(1))]
      partition_index_data <- list(reform_bar2_vals = p_index_data_vals, reform_bar2_xlev = p_index_data_xlev)
      # partition_index_data <- as.vector(model.matrix(reform_bar2, data))
    }
    p_index_newdata_mx <- model.matrix(reform_bar2, model.frame(reform_bar2, newdata, na.action = na.pass, xlev = partition_index_data$reform_bar2_xlev))
    p_index_newdata_names <- colnames(p_index_newdata_mx)
    p_index_newdata_split <- split(p_index_newdata_mx, seq_len(NROW(p_index_newdata_mx)))
    partition_index_newdata <- p_index_newdata_names[vapply(p_index_newdata_split, function(y) which(as.logical(y)), numeric(1))]
    # partition_index_newdata <- as.vector(model.matrix(reform_bar2, newdata)
    partition_index <- vapply(partition_index_newdata, function(x) ifelse(x == partition_index_data$reform_bar2_vals, 1, 0), numeric(length(partition_index_data$reform_bar2_vals)))
    t_partition_index <- Matrix(t(partition_index), sparse = TRUE)
  }
  t_partition_index
}
