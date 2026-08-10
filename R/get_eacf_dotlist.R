#' Get empirical autocovariance dotlist
#'
#' @param ... Additional arguments to \code{eacf()}
#' @param max_halfdist
#'
#' @return An eacf dotlist
#'
#' @noRd
get_eacf_dotlist <- function(..., max_halfdist) {
  # collects user-supplied ... args and fills in bins/cutoff defaults so
  # downstream eacf() calls always receive both, whether or not the caller set them
  # storing dotlist and setting defaults for eacf
  dotlist <- list(...)

  if (!("bins" %in% names(dotlist))) {
    dotlist$bins <- 15
  }

  if (!("cutoff" %in% names(dotlist))) {
    dotlist$cutoff <- max_halfdist
  }

  # make dotlist eacf
  dotlist_eacf <- list(bins = dotlist$bins, cutoff = dotlist$cutoff)
}


#' Fill in default plotting arguments for an empirical autocovariance plot
#'
#' @param x An eacf object (as returned by \code{get_eacf()} or \code{get_eacf_cloud()})
#' @param dotlist A list of user-supplied plotting arguments
#' @param cloud Whether \code{x} is a cloud (unbinned) eacf object
#'
#' @return \code{dotlist} with defaults filled in for any of \code{main},
#'   \code{xlab}, \code{ylab}, \code{pch}, \code{cex}, and \code{ylim} not
#'   already supplied
#'
#' @noRd
get_eacf_dotlist_defaults <- function(x, dotlist, cloud) {
  # fills in base::plot()-style defaults (title, axis labels, point size/shape,
  # y limits) for whichever of these the caller did not already specify
  names_dotlist <- names(dotlist)

  # set defaults
  if (!"main" %in% names_dotlist) {
    dotlist$main <- "Empirical Autocovariance"
    if (cloud) dotlist$main <- paste0(dotlist$main, " (Cloud)")
  }

  if (!"xlab" %in% names_dotlist) {
    dotlist$xlab <- "Distance"
  }

  if (!"ylab" %in% names_dotlist) {
    dotlist$ylab <- "Autocovariance"
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
    # include zero if not in limits

    # na.rm = TRUE throughout: bins beyond the data's actual extent (e.g. a
    # cutoff larger than any observed pairwise distance) have no pairs and so
    # a NA acov: all()/max()/min() must ignore those to still classify sign
    # and find a finite ylim from the real bins
    ## all greater than zero (positive)
    if (all(x$acov > 0, na.rm = TRUE)) {
      dotlist$ylim <- c(0, 1.1 * max(x$acov, na.rm = TRUE))
    }

    ## all less than zero (negative)
    if (all(x$acov < 0, na.rm = TRUE)) {
      dotlist$ylim <- c(1.1 * min(x$acov, na.rm = TRUE), 0)
    }
  }

  dotlist
}
