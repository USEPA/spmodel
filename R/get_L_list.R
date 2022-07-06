#' Get relevant L lists for anova
#'
#' @param assign_index A single assign value from the model matrix
#' @param assign_indices The assign values from the model matrix
#'
#' @return L lists for anova
#'
#' @noRd
get_L_list <- function(assign_index, assign_indices) {
  assign_vals <- which(assign_indices == assign_index)
  L_vectors <- lapply(assign_vals, get_L_vector, assign_indices)
  do.call(rbind, L_vectors)
}

get_L_vector <- function(assign_val, assign_indices) {
  L_vector <- matrix(0, nrow = 1, ncol = length(assign_indices))
  L_vector[, assign_val] <- 1
  L_vector
}
