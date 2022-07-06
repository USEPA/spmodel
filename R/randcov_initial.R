#' Create a random effects covariance parameter initial object
#'
#' @description Create a random effects (co)variance parameter initial object that specifies
#'   initial and/or known values to use while estimating random effect variances
#'   with [splm()] or [spautor()].
#'
#' @param ... Arguments to \code{randcov_params()}.
#' @param known A character vector indicating which random effect variances are to be
#'   assumed known. The value \code{"given"} is shorthand for assuming all
#'   random effect variances given to \code{randcov_initial()} are assumed known.
#'
#' @details A random effect is specified as \eqn{Zu}, where \eqn{Z} is the random
#'   effects design matrix and \code{u} is the random effect. The covariance of
#'   \eqn{Zu} is \eqn{\sigma 2 ZZ^T}, where \eqn{\sigma 2} is the random effect
#'   variance, and \eqn{Z^T} is the transpose of \eqn{Z}.
#'
#' @return A list with two elements: \code{initial} and \code{is_known}.
#'   \code{initial} is a named numeric vector indicating the random effect variances
#'   with specified initial and/or known values. \code{is_known} is a named
#'   logical vector indicating whether the random effect variances in
#'   \code{initial} are known or not.
#'
#' @export
#'
#' @examples
#' print(randcov_initial(group = 1))
#' print(randcov_initial(group = 1, known = "group"))
randcov_initial <- function(..., known) {
  randcov_params_given <- randcov_params(...)
  if (missing(known)) {
    is_known <- rep(FALSE, length(randcov_params_given))
  } else {
    if (identical(known, "given")) {
      is_known <- rep(TRUE, length(randcov_params_given))
    } else {
      is_known <- names(randcov_params_given) %in% known
    }
  }
  names(is_known) <- names(randcov_params_given)

  # error if NA and known
  randcov_NA <- which(is.na(randcov_params_given))
  if (any(is_known[randcov_NA])) {
    stop("randcov_initial values cannot be NA and known.", call. = FALSE)
  }

  new_randcov_initial <- list(initial = randcov_params_given, is_known = is_known)
  new_randcov_initial
}
