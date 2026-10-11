# Earlier versions averaged nonzero partition BLUPs into one value per level.
# That summarizes distinct effects but is not a partition-level BLUP, which is now returned.
get_fitted_randcov <- function(data_object, randcov_params, SigInv_r_list) {
  if (is.null(names(randcov_params))) return(NULL)
  partition <- if ("partition_factor_initial" %in% names(data_object)) data_object$partition_factor_initial else data_object$partition_factor
  design <- get_randcov_design(data_object$obdata, randcov_params, partition)
  index <- if (is.null(data_object$local_index)) list(seq_len(NROW(data_object$obdata))) else
    unname(split(seq_len(NROW(data_object$obdata)), data_object$local_index))
  total <- numeric(length(design$variance))
  count <- integer(length(total))
  for (k in seq_along(index)) {
    value <- design$variance * as.numeric(crossprod(design$observed[index[[k]], , drop = FALSE], SigInv_r_list[[k]]))
    total <- total + value
    count <- count + (value != 0)
  }
  value <- total / pmax(count, 1L)
  out <- lapply(names(randcov_params), function(term) {
    rows <- which(design$info$term == term)
    ans <- setNames(value[rows], design$info$level[rows])
    if (is.null(partition)) {
      # Preserve the existing zero entries for declared but unobserved levels.
      Z <- if (!is.null(data_object$randcov_list)) data_object$randcov_list[[1L]][[term]]$Z else
        get_randcov_Z(term, data_object$obdata, ZZt = FALSE)$Z
      full <- setNames(numeric(NCOL(Z)), colnames(Z))
      full[names(ans)] <- ans
      ans <- full
    }
    ans
  })
  setNames(out, names(randcov_params))
}
