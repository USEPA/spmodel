#' Transform coordinates to accommodate anisotropy
#'
#' @param rotate The rotate value (between 0 and pi radians)
#' @param scale The scale value for the minor axis (between 0 and 1)
#' @param xcoord_val The x-coordinate value (Euclidean isotropic)
#' @param ycoord_val The y-coordinate value (Euclidean isotropic)
#'
#' @return New coordinates
#'
#' @noRd
transform_anis <- function(data, xcoord, ycoord, rotate, scale) {
  rotate_clockwise <- matrix(c(cos(rotate), sin(rotate), -sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
  scale_yaxis <- matrix(c(1, 0, 0, 1 / scale), nrow = 2, ncol = 2, byrow = TRUE)
  coords <- rbind(data[[xcoord]], data[[ycoord]])
  new_coords <- (scale_yaxis %*% rotate_clockwise) %*% coords
  list(xcoord_val = new_coords[1, ], ycoord_val = new_coords[2, ])
}

# transform_anis <- function(data, xcoord, ycoord, xcoord_val, ycoord_val, rotate, scale) {
#   rotate_clockwise <- matrix(c(cos(rotate), sin(rotate), -sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
#   scale_yaxis <- matrix(c(1, 0, 0, 1 / scale), nrow = 2, ncol = 2, byrow = TRUE)
#   if (missing(xcoord) && missing(ycoord)) {
#     coords <- rbind(xcoord_val, ycoord_val)
#   } else {
#     coords <- rbind(data[[xcoord]], data[[ycoord]])
#   }
#   new_coords <- (scale_yaxis %*% rotate_clockwise) %*% coords
#   list(xcoord_val = new_coords[1, ], ycoord_val = new_coords[2, ])
# }
