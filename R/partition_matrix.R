#' Create a partition matrix
#'
#' @param partition_factor A partition factor (formula)
#' @param data Data
#'
#' @return A partition matrix
#'
#' @noRd
partition_matrix <- function(partition_factor = NULL, data) {
  # a partition factor restricts spatial dependence to within groups: the
  # returned matrix is 1 for pairs of observations in the same group and 0
  # otherwise, and gets multiplied elementwise into the covariance matrix
  # later to zero out covariance between different groups
  if (is.null(partition_factor)) {
    partition_matrix_val <- NULL
  } else {
    # finding the formula
    partition_formula <- reformulate(labels(terms(partition_factor)), intercept = FALSE)
    # use regular contrasts here so matrix all zeros and ones
    partition_model_frame <- model.frame(partition_formula, data)
    if (length(unique(as.character(unlist(partition_model_frame)))) == 1) {
      partition_model_val <- Matrix::Matrix(matrix(1, nrow = NROW(partition_model_frame), ncol = 1), sparse = TRUE)
    } else {
      # built directly as a sparse indicator from the group-label factor
      # (the interaction of the partition variables) rather than via a dense
      # model.matrix(), which would allocate an n x nlevels dense matrix
      # before ever going sparse -- expensive when there are many levels
      group_label <- interaction(partition_model_frame, drop = FALSE)
      partition_model_val <- Matrix::t(Matrix::fac2sparse(group_label, drop.unused.levels = FALSE))
    }
    # tcrossprod of a 0/1 group-indicator matrix with itself gives, at
    # position (i, j), the dot product of observation i's and j's indicator
    # rows -- 1 if they share a group (both have a 1 in the same column) and
    # 0 otherwise, which is exactly the desired n x n group-membership matrix
    partition_matrix_val <- tcrossprod(partition_model_val, partition_model_val)
  }
  partition_matrix_val
}

#' Get each observation's partition factor group label
#'
#' Big data efficient companion to \code{partition_matrix()}: returns the length-n
#' vector of group labels (the interaction of the partition variables)
#' instead of building the full n x n group-membership matrix.
#'
#' @param partition_factor A partition factor (formula)
#' @param data Data
#'
#' @return A factor of group labels, one per row of \code{data}, or
#'   \code{NULL} if \code{partition_factor} is \code{NULL}.
#'
#' @noRd
partition_group <- function(partition_factor = NULL, data) {
  if (is.null(partition_factor)) {
    return(NULL)
  }
  partition_formula <- reformulate(labels(terms(partition_factor)), intercept = FALSE)
  partition_model_frame <- model.frame(partition_formula, data)
  interaction(partition_model_frame, drop = FALSE)
}
