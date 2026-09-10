#' Restore the original y-coordinate name/values for 1-dimensional coordinates
#'
#' @param data_object The data object
#'
#' @return \code{data_object} with its (internally-generated) y-coordinate
#'   column renamed and repopulated with the user's original y-coordinate
#'   name/values, when the data are 1-dimensional (\code{dim_coords == 1})
#'
#' @noRd
replace_data_object_dimcoords1 <- function(data_object) {
  # replace original coordinate
  # for 1-dimensional data, internal computations still need a y-coordinate,
  # so one was generated (and obdata/newdata were given a placeholder ycoord
  # column) upstream; ycoord_orig_name/_val being non-NULL signals that
  # substitution happened and needs to be undone before returning to the user
  if (!is.null(data_object$ycoord_orig_name)) {
    # rename the placeholder column back to the user's original y-coordinate name
    names(data_object$obdata)[which(names(data_object$obdata) == data_object$ycoord)] <- as.character(data_object$ycoord_orig_name)
    # and restore its original values, indexed to the observed (non-missing) rows
    data_object$obdata[[data_object$ycoord_orig_name]] <- data_object$ycoord_orig_val[data_object$observed_index]
    if (!is.null(data_object$newdata)) {
      names(data_object$newdata)[which(names(data_object$newdata) == data_object$ycoord)] <- as.character(data_object$ycoord_orig_name)
      # newdata rows correspond to the missing_index subset instead
      data_object$newdata[[data_object$ycoord_orig_name]] <- data_object$ycoord_orig_val[data_object$missing_index]
    }
    data_object$ycoord <- data_object$ycoord_orig_name
  }
  data_object
}
