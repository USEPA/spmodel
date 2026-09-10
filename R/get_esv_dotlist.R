#' Get empirical semivariogram dotlist
#'
#' @param ... Additional arguments to \code{esv()}
#' @param max_halfdist
#'
#' @return An esv dotlist
#'
#' @noRd
get_esv_dotlist <- function(..., max_halfdist) {
  # collects user-supplied ... args and fills in bins/cutoff/robust defaults so
  # downstream esv() calls always receive all three, whether or not the caller set them
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


#' Fill in default plotting arguments for an empirical semivariogram plot
#'
#' @param x An esv object (as returned by \code{get_esv()}, \code{get_esv_robust()}, or \code{get_esv_cloud()})
#' @param dotlist A list of user-supplied plotting arguments
#' @param cloud Whether \code{x} is a cloud (unbinned) esv object
#'
#' @return \code{dotlist} with defaults filled in for any of \code{main},
#'   \code{xlab}, \code{ylab}, \code{pch}, \code{cex}, and \code{ylim} not
#'   already supplied
#'
#' @noRd
get_esv_dotlist_defaults <- function(x, dotlist, cloud) {
  # fills in base::plot()-style defaults (title, axis labels, point size/shape,
  # y limits) for whichever of these the caller did not already specify
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
    # scale point size by number of pairs (np) in each bin, rescaled to [1, 3],
    # so bins backed by more pairs (more reliable estimates) plot larger
    dotlist$cex <- (x$np - min(x$np)) / (max(x$np) - min(x$np)) * 2 + 1
  }

  if (!"ylim" %in% names_dotlist) {
    # na.rm = TRUE: bins beyond the data's actual extent (e.g. a cutoff larger
    # than any observed pairwise distance) have no pairs and so a NA gamma;
    # max() must ignore those to still find a finite ylim from the real bins
    dotlist$ylim <- c(0, 1.1 * max(x$gamma, na.rm = TRUE))
  }

  dotlist
}
