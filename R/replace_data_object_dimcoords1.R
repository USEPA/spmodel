replace_data_object_dimcoords1 <- function(data_object) {
  # replace original coordinate
  if (!is.null(data_object$ycoord_orig_name)) {
    names(data_object$obdata)[which(names(data_object$obdata) == data_object$ycoord)] <- as.character(data_object$ycoord_orig_name)
    data_object$obdata[[data_object$ycoord_orig_name]] <- data_object$ycoord_orig_val[data_object$observed_index]
    if (!is.null(data_object$newdata)) {
      names(data_object$newdata)[which(names(data_object$newdata) == data_object$ycoord)] <- as.character(data_object$ycoord_orig_name)
      data_object$newdata[[data_object$ycoord_orig_name]] <- data_object$ycoord_orig_val[data_object$missing_index]
    }
    data_object$ycoord <- data_object$ycoord_orig_name
  }
  data_object
}
