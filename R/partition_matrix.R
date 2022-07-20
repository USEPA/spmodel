#' Create a partition matrix
#'
#' @param partition_factor A partition factor (formula)
#' @param data Data
#'
#' @return A partition matrix
#'
#' @noRd
partition_matrix <- function(partition_factor = NULL, data) {
  if (is.null(partition_factor)) {
    partition_matrix_val <- NULL
  } else {
    # finding the formula
    partition_formula <- reformulate(labels(terms(partition_factor)), intercept = FALSE)
    # use regular contrasts here so matrix all zeros and ones
    partition_model_frame <- model.frame(partition_formula, data)
    if (length(unique(as.character(unlist(partition_model_frame)))) == 1) {
      partition_model_val <- Matrix::Matrix(as.matrix(partition_model_frame), sparse = TRUE)
    } else {
      partition_model_val <- Matrix::Matrix(model.matrix(partition_formula, data), sparse = TRUE)
    }
    partition_matrix_val <- tcrossprod(partition_model_val, partition_model_val)
  }
  partition_matrix_val
}
