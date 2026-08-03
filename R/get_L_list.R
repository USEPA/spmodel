#' Get relevant L lists for anova
#'
#' @param assign_index A single assign value from the model matrix
#' @param assign_indices The assign values from the model matrix
#'
#' @return L lists for anova
#'
#' @noRd
get_L_list <- function(assign_index, assign_indices) {
  # assign_indices maps each model matrix column to the model term it belongs
  # to (as produced by attr(model.matrix(...), "assign")); find every column
  # belonging to this particular term
  assign_vals <- which(assign_indices == assign_index)
  # build one indicator row per column, then stack them into the L matrix
  # used to test this term's joint contrast (L %*% beta = 0) in anova()
  L_vectors <- lapply(assign_vals, get_L_vector, assign_indices)
  do.call(rbind, L_vectors)
}

#' Get a single indicator row vector for anova's L matrix
#'
#' @param assign_val The column index (into the model matrix) to indicate
#' @param assign_indices The assign values from the model matrix
#'
#' @return A single-row matrix of zeros with a one in column \code{assign_val}
#'
#' @noRd
get_L_vector <- function(assign_val, assign_indices) {
  L_vector <- matrix(0, nrow = 1, ncol = length(assign_indices))
  # a single 1 at position assign_val picks out that coefficient when
  # multiplied against the full coefficient vector (L %*% beta)
  L_vector[, assign_val] <- 1
  L_vector
}
