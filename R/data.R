#' Data on heavy metals in mosses near a mining road in Alaska, USA
#'
#' @description Data on heavy metals in mosses near a mining road in Alaska, USA.
#'
#' @format An \code{sf} object with 365 rows and 10 columns:
#'
#' \itemize{
#'   \item{sample: }{A factor with a sample identifier. Some samples were
#'     replicated in the field or laboratory. As a result, there are 318 unique
#'     sample identifiers.}
#'   \item{field_rep: }{A factor representing field replicate. Takes values \code{1}
#'     and \code{2}.}
#'   \item{lab_rep: }{A factor representing laboratory replicate. Takes values \code{1}
#'     and \code{2}.}
#'   \item{year: }{A factor representing year. Takes values \code{2001} and \code{2006}.}
#'   \item{sideroad: }{A factor representing direction relative to the haul road.
#'     Takes values \code{N} (north of the haul road) and \code{S} (south
#'     of the haul road).}
#'   \item{log_dist2road: }{The log of distance (in meters) to the haul road.}
#'   \item{log_Zn: }{The log of zinc concentration in moss tissue (mg/kg).}
#'   \item{geometry: }{\code{POINT} geometry representing coordinates in an Alaksa
#'     Albers projection (EPSG: 3338).}
#' }
#' @source Data were obtained from Peter Neitlich and Linda Hasselbach of the National
#'   Park Service.  Data were used in the publications listed in References.
#' @references
#' Neitlich, P.N., Ver Hoef, J.M., Berryman, S. D., Mines, A., Geiser, L.H.,
#'   Hasselbach, L.M., and Shiel, A. E. 2017. Trends in Spatial Patterns of Heavy
#'   Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul
#'   Road, Alaska, 2001-2006. PLOS ONE 12(5):e0177936 DOI:10.1371/journal.pone.0177936
#'   (\href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177936}{links
#'   to: Trends in Spatial Patterns of Heavy Metal Deposition on National
#'   Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001-2006.})
#'
#' Hasselbach, L., Ver Hoef, J.M., Ford, J., Neitlich, P., Berryman, S., Wolk B.
#'   and Bohle, T. 2005. Spatial Patterns of Cadmium, Lead and Zinc Deposition
#'   on National Park Service Lands in the Vicinity of Red Dog Mine, Alaska.
#'   Science of the Total Environment 348: 211-230.
#'   (\href{https://www.sciencedirect.com/science/article/pii/S0048969705000082}{links
#'   to: Spatial Patterns of Cadmium, Lead and Zinc Deposition on National
#'   Park Service Lands in the Vicinity of Red Dog Mine, Alaska.})
"moss"



#' Data for a caribou forage experiment
#'
#' @description Data for a caribou forage experiment.
#'
#' @format A \code{tibble} with 30 rows and 5 columns:
#' \itemize{
#'   \item{water: }{A factor representing whether water was added. Takes values
#'     \code{N} (no water added) and \code{Y} (water added).}
#'   \item{tarp: }{A factor representing tarp cover. Takes values \code{clear}
#'     (a clear tarp), \code{shade} (a shade tarp), and \code{none} (no tarp).}
#'   \item{z: }{The percentage of nitrogen.}
#'   \item{x: }{The x-coordinate.}
#'   \item{y: }{The y-coordinate.}
#' }
#' @source These data were provided by Elizabeth Lenart of the Alaska Department
#'   of Fish and Game.  The data were used in the publication listed in References.
#' @references Lenart, E.A., Bowyer, R.T., Ver Hoef, J.M. and Ruess, R.W. 2002.
#'   Climate Change and Caribou: Effects of Summer Weather on Forage. Canadian
#'   Journal of Zoology 80: 664-678.
#'   (\href{https://www.nrcresearchpress.com/doi/abs/10.1139/z02-034#.XvweB3VKjmF}{links
#'   to: Climate Change and Caribou: Effects of Summer Weather on Forage})
"caribou"

#' Data on estimated harbor seal trends in southeast Alaska, USA
#'
#' @description Data on estimated harbor seal abundance trends in southeast Alaska, USA.
#'
#' @format A \code{sf} object with 62 rows and 2 columns:
#' \describe{
#'   \item{log_trend: }{The log of the estimated harbor seal abundance trends.}
#'   \item{geometry: }{\code{POLYGON} geometry representing polygons in an Alaska
#'     Albers projection (EPSG: 3338).}
#' }
#' @source These data were collected by the Polar Ecosystem Program of the Marine
#'   Mammal Laboratory of the Alaska Fisheries Science Center of NOAA Fisheries.
#'   The data were used in the publication listed in References.
#' @references
#' Ver Hoef, J.M., Peterson, E. E., Hooten, M. B., Hanks, E. M., and Fortin, M.-J. 2018.
#'   Spatial Autoregressive Models for Statistical Inference from Ecological Data.
#'   Ecological Monographs, 88: 36-59. DOI: 10.1002/ecm.1283.
#'   \href{https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1283}{links
#'   to: Spatial Autoregressive Models for Statistical Inference from Ecological Data.}.
"seal"

#' Data on sulfate atmospheric deposition in the conterminous USA
#'
#' @description Data on sulfate atmospheric deposition in the conterminous USA.
#'
#' @format An \code{sf} object with 197 rows and 2 columns.
#' \describe{
#'   \item{sulfate: }{Total wet deposition sulfate in kilograms per hectare.}
#'   \item{geometry: }{\code{POINT} geometry representing coordinates in a
#'     Conus Albers projection (EPSG: 5070).}
#' }
#' @source
#' These data were used in the publication listed in References. Data were downloaded from the website
#' (\href{http://nadp.slh.wisc.edu/NTN/}{links to: National Atmospheric Deposition Program}).
#' @references
#' Zimmerman, D.L. (1994). Statistical analysis of spatial data. Pages 375-402 in
#'   \emph{ Statistical Methods for Physical Science}, J. Stanford and
#'   S. Vardeman (eds.), Academic Press: New York.
#'   \href{https://www.elsevier.com/books/statistical-methods-for-physical-science/stanford/978-0-12-475973-2}{links
#'   to: Statistical analysis of spatial data}.
"sulfate"

#' Prediction locations for sulfate atmospheric deposition in the conterminous USA
#'
#' @description Prediction locations for sulfate atmospheric deposition in the conterminous USA.
#'
#' @format An \code{sf} object with 197 rows and 1 column.
#' \describe{
#'   \item{geometry: }{\code{POINT} geometry representing coordinates in a
#'     Conus Albers projection (EPSG: 5070).}
#' }
#' @source
#' These data were used in the publication listed in References. Data were downloaded from the website
#' (\href{http://nadp.slh.wisc.edu/NTN/}{links to: National Atmospheric Deposition Program}).
#' @references
#' Zimmerman, D.L. (1994). Statistical analysis of spatial data. Pages 375-402 in
#'   \emph{ Statistical Methods for Physical Science}, J. Stanford and
#'   S. Vardeman (eds.), Academic Press: New York.
#'   \href{https://www.elsevier.com/books/statistical-methods-for-physical-science/stanford/978-0-12-475973-2}{links
#'   to: Statistical analysis of spatial data}.
"sulfate_preds"
