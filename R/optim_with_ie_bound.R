#' Call optim with an optional independent-error lower bound
#'
#' @noRd
optim_with_ie_bound <- function(optim_args, optim_dotlist, lower = NULL) {
  if (is.null(lower)) {
    return(do.call("optim", c(optim_args, optim_dotlist)))
  }
  parameter <- names(optim_args$par)
  upper <- optim_dotlist$upper
  if (length(parameter) != 1L || !parameter %in% c("ie_log", "ie_prop_logodds") ||
    length(lower) != 1L || !is.finite(lower) || length(upper) != 1L || !is.finite(upper) || lower >= upper) {
    stop("Invalid bounds for the independent-error variance search.", call. = FALSE)
  }
  optim_dotlist$lower <- lower
  do.call("optim", c(optim_args, optim_dotlist))
}

#' Lower search bound for a sole unknown independent-error variance
#'
#' @noRd
ie_only_lower <- function(par, spcov_initial, diagtol, profiled = FALSE) {
  parameter <- if (profiled) "ie_prop_logodds" else "ie_log"
  if (!identical(names(par), parameter) || inherits(spcov_initial, c("none", "car", "sar"))) {
    return(NULL)
  }
  # The profiled plateau is scale-invariant only with zero absolute tolerance.
  if (profiled && diagtol != 0) return(NULL)
  spcov <- spcov_initial$initial
  spcov[["ie"]] <- 0
  if (profiled) spcov[["de"]] <- 1
  floor <- spcov_ie_stabilized(spcov, diagtol)
  if (isTRUE(floor == 0)) return(NULL)
  log(floor)
}
