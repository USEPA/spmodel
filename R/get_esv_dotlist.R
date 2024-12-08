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

  if (!("robust" %in% names(dotlist))) {
    dotlist$robust <- FALSE
  }

  # make dotlist esv
  dotlist_esv <- list(robust = dotlist$robust, bins = dotlist$bins, cutoff = dotlist$cutoff)
}


get_esv_dotlist_defaults <- function(x, dotlist, cloud) {

  names_dotlist <- names(dotlist)

  # set defaults
  if (!"main" %in% names_dotlist) {
    dotlist$main <- "Empirical Semivariogram"
    if (cloud) dotlist$main <- paste0(dotlist$main, " (Cloud)")
  }

  if (!"xlab" %in% names_dotlist) {
    dotlist$xlab <- "Distance"
  }

  if (!"ylab" %in% names_dotlist) {
    dotlist$ylab <- "Semivariance"
  }

  if (!cloud && !"pch" %in% names_dotlist) {
    dotlist$pch <- 19
  }

  if (!cloud && !"cex" %in% names_dotlist) {
    dotlist$cex <- (x$np - min(x$np)) / (max(x$np) - min(x$np)) * 2 + 1
  }

  if (!"ylim" %in% names_dotlist) {
    dotlist$ylim <- c(0, 1.1 * max(x$gamma))
  }

  dotlist
}
