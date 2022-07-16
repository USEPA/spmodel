#' Turn an \code{sf} object into a data frame
#'
#' @param data An \code{sf} object with all \code{POINT} geometries or geometries
#'   coercible to \code{POINT} geometries via sf::st_centriod().
#'
#' @return A data frame with coordinates added as columns
#'
#' @noRd
sf_to_df <- function(data) {
  ### find the sf geometry column location
  geometry_column <- which(colnames(data) == attr(data, "sf_column"))
  ### see if the geometry column has sfc_point class (a point geometry)
  if (inherits(data[[geometry_column]], "sfc_POINT")) {
    # save coordinates
    coords <- sf::st_coordinates(data)
    # add coordinates to data
    data$xcoord <- coords[, "X"]
    data$ycoord <- coords[, "Y"]
    # drop geometry from data
    data <- sf::st_drop_geometry(data)
  } else {
    ### return an error if all geometries are not POINT (user must use
    ### sf::st_cast() or an equivalent)
    stop("All sf geometries must be POINT geometries or able to be coerced to POINT geometries by sf::st_centroid()")
  }
  data
}

# sf_to_df <- function(data) {
#   ### find the sf geometry column location
#   geometry_column <- which(colnames(data) == attr(data, "sf_column"))
#   ### see if the geometry column has sfc_point class (a point geometry)
#   if (inherits(data[[geometry_column]], "sfc_POINT")) {
#     ### save the coordinates
#     #### (structured as xcoord1, ycoord1, xcoord2, ycoord2, ... , ycoordn)
#     coords <- unlist(data[[geometry_column]])
#     ### redefine data as a data frame
#     data <- as.data.frame(data)
#     ### drop the geometry column
#     data <- data[, -geometry_column, drop = FALSE]
#     ### take the xcoordinate from coords and name it xcoord
#     data$xcoord <- coords[seq(from = 1, to = length(coords), by = 2)]
#     ### take the ycoordinate from coords and name it ycoord
#     data$ycoord <- coords[seq(from = 2, to = length(coords), by = 2)]
#   } else {
#     ### return an error if all geometries are not POINT (user must use
#     ### sf::st_cast() or an equivalent)
#     stop("All geometries in the sf object must be POINT, POLYGON, or MULTIPOLYGON geometries")
#   }
#   data
# }

# # Turn an \code{sp} object into a data frame
# #
# # @param data A \code{SpatialPointsDataFrame} \code{sp} object
# #
# # @return A data frame with coordinates added as columns
# #
# # @noRd
# sp_to_df <- function(data) {
#   ## convert sp to data frame (point geometry)
#   ### see if data has SpatialPointsDataFrame class (a point geometry)
#   ### take coordinate slot in data
#   if (inherits(data, "SpatialPointsDataFrame")) {
#     coords <- data@coords
#     ### take data slot in data and redefine data (as a data frame)
#     data <- data@data
#     ### take the xcoordinate from coords and name it xcoord and ycoord
#     data$xcoord <- coords[, 1]
#     data$ycoord <- coords[, 2]
#     return(data)
#   } else {
#     stop("sp object must have class SpatialPointsDataFrame")
#   }
# }
