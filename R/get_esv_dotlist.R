#' Get empirical semivariogram dotlist
#'
#' @param ... Additional arguments to \code{esv()}
#' @param max_halfdist
#'
#' @return An esv dotlist
#'
#' @noRd
get_esv_dotlist <- function(..., max_halfdist) {
  # storing dotlist and setting defaults for esv
  dotlist <- list(...)

  if (!("bins" %in% names(dotlist))) {
    dotlist$bins <- 15
  }

  if (!("cutoff" %in% names(dotlist))) {
    dotlist$cutoff <- max_halfdist
  }
  # make dotlist esv
  dotlist_esv <- list(bins = dotlist$bins, cutoff = dotlist$cutoff)
}
