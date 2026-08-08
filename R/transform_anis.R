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
# geometric anisotropy is modeled by rotating coordinates so the spatial
# process's major axis of dependence aligns with the x-axis, then stretching
# the (now axis-aligned) minor axis so an anisotropic process becomes
# isotropic in the transformed space, where ordinary Euclidean distance can
# be used again
transform_anis <- function(data, xcoord, ycoord, rotate, scale) {
  # clockwise rotation matrix by "rotate" radians
  rotate_clockwise <- matrix(c(cos(rotate), sin(rotate), -sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
  # 1/scale on the y-axis stretches the minor axis back to match the major
  # axis's range dependence (scale is between 0 and 1)
  scale_yaxis <- matrix(c(1, 0, 0, 1 / scale), nrow = 2, ncol = 2, byrow = TRUE)
  coords <- rbind(data[[xcoord]], data[[ycoord]])
  # rotate first, then scale, applied as a single combined linear map
  new_coords <- (scale_yaxis %*% rotate_clockwise) %*% coords
  list(xcoord_val = new_coords[1, ], ycoord_val = new_coords[2, ])
}

#' The inverse of transforming coordinates to accommodate anisotropy
#'
#' @return new coordinates
#'
#' @noRd
# undoes transform_anis(): unscale then rotate counterclockwise, i.e. applies
# the inverse operations in reverse order, to map isotropic coordinates back
# to their original anisotropic space
transform_anis_inv <- function(data, xcoord, ycoord, rotate, scale) {
  rotate_clockwise_inv <- matrix(c(cos(rotate), -sin(rotate), sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
  scale_yaxis_inv <- matrix(c(1, 0, 0, scale), nrow = 2, ncol = 2, byrow = TRUE)
  coords <- rbind(data[[xcoord]], data[[ycoord]])
  new_coords <- (rotate_clockwise_inv %*% scale_yaxis_inv) %*% coords
  list(xcoord_val = new_coords[1, ], ycoord_val = new_coords[2, ])
}

transform_anis2 <- function(xcoord_val, ycoord_val, rotate, scale) {
  rotate_clockwise <- matrix(c(cos(rotate), sin(rotate), -sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
  scale_yaxis <- matrix(c(1, 0, 0, 1 / scale), nrow = 2, ncol = 2, byrow = TRUE)
  coords <- rbind(xcoord_val, ycoord_val)
  new_coords <- (scale_yaxis %*% rotate_clockwise) %*% coords
  list(xcoord_val = new_coords[1, ], ycoord_val = new_coords[2, ])
}
