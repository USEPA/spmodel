#' Turn an \code{sf} object into a data frame
#'
#' @param data An \code{sf} object with all \code{POINT} geometries or geometries
#'   coercible to \code{POINT} geometries via sf::st_centroid().
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
    data$.xcoord <- coords[, "X"]
    data$.ycoord <- coords[, "Y"]
    # drop geometry from data
    data <- sf::st_drop_geometry(data)
  } else {
    ### return an error if all geometries are not POINT (user must use
    ### sf::st_cast() or an equivalent)
    stop("All sf geometries must be POINT geometries or able to be coerced to POINT geometries by sf::st_centroid()")
  }
  data
}
