#' Get empirical autocovariance dotlist
#'
#' @param ... Additional arguments to \code{eacf()}
#' @param max_halfdist
#'
#' @return An eacf dotlist
#'
#' @noRd
get_eacf_dotlist <- function(..., max_halfdist) {
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


get_eacf_dotlist_defaults <- function(x, dotlist, cloud) {

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
    dotlist$cex <- (x$np - min(x$np)) / (max(x$np) - min(x$np)) * 2 + 1
  }

  if (!"ylim" %in% names_dotlist) {

    # include zero if not in limits

    ## all greater than zero (positive)
    if (all(x$acov > 0)) {
      dotlist$ylim <- c(0, 1.1 * max(x$acov))
    }

    ## all less than zero (negative)
    if (all(x$acov < 0)) {
      dotlist$ylim <- c(1.1 * min(x$acov), 0)
    }
  }

  dotlist
}
