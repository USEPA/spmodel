#' @export
# use generics to export tidy
generics::tidy

#' @export
# use generics to export glance
generics::glance

#' @export
# use generics to export augment
generics::augment

#' Compute \eqn{diag(\mathbf{X}\mathbf{V}\mathbf{X}^\top)} without forming the
#' full product
#'
#' @param X An \eqn{n \times p} matrix
#' @param V A \eqn{p \times p} matrix
#'
#' @details Only the diagonal is ever wanted, and the \eqn{i}th diagonal entry
#'   is the scalar \eqn{\mathbf{x}_i^\top \mathbf{V} \mathbf{x}_i}. Since
#'   \eqn{(\mathbf{X}\mathbf{V})} has \eqn{\mathbf{x}_i^\top \mathbf{V}} as its
#'   \eqn{i}th row, multiplying it elementwise by \eqn{\mathbf{X}} and summing
#'   along each row recovers that scalar. The result is identical to
#'   \code{diag(X \%*\% V \%*\% t(X))}, but costs \eqn{O(np^2)} operations and
#'   \eqn{O(np)} memory rather than \eqn{O(n^2p)} and \eqn{O(n^2)}.
#'
#' @return A numeric vector of length \code{NROW(X)}
#'
#' @noRd
get_diag_XVXt <- function(X, V) {
  as.numeric(rowSums((X %*% V) * X))
}

#' Standard error of the fitted mean at the observed locations
#'
#' @param object A fitted model object
#'
#' @details This is \eqn{\sqrt{diag(\mathbf{X}(\mathbf{X}^\top
#'   \boldsymbol{\Sigma}^{-1}\mathbf{X})^{-1}\mathbf{X}^\top)}}, the standard
#'   error attaching to \eqn{\mathbf{X}\hat{\boldsymbol{\beta}}}, since
#'   \code{vcov()} already returns \eqn{(\mathbf{X}^\top
#'   \boldsymbol{\Sigma}^{-1}\mathbf{X})^{-1}}. It is a confidence (fitted mean)
#'   standard error, not a prediction standard error, and matches what
#'   \code{predict(interval = "confidence")} reports. For generalized linear
#'   models it is on the link scale, again matching \code{predict()}.
#'
#' @return A vector of standard errors, one per observed location
#'
#' @noRd
get_se_fitted_mean <- function(object) {
  sqrt(get_diag_XVXt(model.matrix(object), vcov(object)))
}

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
