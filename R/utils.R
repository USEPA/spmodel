#' @export
# use generics to export tidy
generics::tidy

#' @export
# use generics to export glance
generics::glance

#' @export
# use generics to export augment
generics::augment

#' Compute the logit (log-odds) transform
#'
#' @param x A value between zero and one
#'
#' @return \code{log(x / (1 - x))}
#'
#' @noRd
# used to map covariance parameters that must stay between zero and one
# (e.g., proportions) onto an unconstrained scale that optim() can search freely
logit <- function(x) {
  if (x < 0 | x > 1) {
    stop("logit argument must be between zero and one", call. = FALSE)
  }
  log(x / (1 - x))
}

#' Compute the expit (inverse logit) transform
#'
#' @param x A value on the log-odds scale
#'
#' @return \code{1 / (1 + exp(-x))}
#'
#' @noRd
# inverse of logit(): converts an optimizer's unconstrained parameter value
# back to the original zero-to-one scale
expit <- function(x) {
  1 / (1 + exp(-x))
}

#' Checklist shown by \code{devtools::release()} before submitting to CRAN
#'
#' @return A character vector of release-checklist questions
#'
#' @noRd
release_questions <- function() {
  c(
    "Have you changed version numbers DESCRIPTION and NEWS?",
    "Are there any loose browser() statements?"
  )
}
