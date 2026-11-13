#' Reconstruct the random formula for specified variance components
#'
#' @param randcov_params Named random effect variances.
#' @param env The generator's evaluation environment.
#' @return A random effects formula with no implicit slope intercepts.
#' @noRd
get_simulation_random <- function(randcov_params, env) {
  if (is.null(names(randcov_params)) || anyNA(names(randcov_params)) ||
      any(!nzchar(trimws(names(randcov_params))))) {
    stop("randcov_params must have a name for each variance component.", call. = FALSE)
  }
  random_terms <- vapply(names(randcov_params), function(label) {
    label <- labels(terms(reformulate(label)))
    if (length(label) != 1L) {
      stop("Each random effect must have only one variance parameter.", call. = FALSE)
    }
    parts <- trimws(strsplit(label, "|", fixed = TRUE)[[1L]])
    if (length(parts) == 1L) {
      val <- paste0("(1 | ", parts[1L], ")")
    } else if (length(parts) == 2L) {
      val <- paste0("(0 + ", parts[1L], " | ", parts[2L], ")")
    } else {
      stop("Invalid random effect variance component name.", call. = FALSE)
    }
    return(val)
  }, character(1))
  init <- randcov_initial(randcov_params, known = "given")
  names(init$initial) <- names(init$is_known) <- random_terms
  validate_randcov_initial(init)
  reformulate(random_terms, env = env)
}
