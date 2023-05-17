#' Create a dispersion parameter initial object
#'
#' @description Create a dispersion parameter initial object that specifies
#'   initial and/or known values to use while estimating the dispersion parameter
#'   with [spglm()] or [spgautor()].
#'
#' @param family The generalized linear model family describing the distribution
#'   of the response variable to be used. \code{"poisson"}, \code{"nbinomial"}, \code{"binomial"},
#'   \code{"beta"}, \code{"Gamma"}, and \code{"inverse.gaussian"}.
#' @param dispersion The value of the dispersion parameter.
#' @param known A character vector indicating whether the dispersion parameter is to be
#'   assumed known. The value \code{"dispersion"} or \code{"given"} is assumes
#'   the dispersion parameter is known.
#'
#' @details The \code{dispersion_initial} list is later passed to [spglm()] or [spgautor()].
#'
#'   The variance function of an individual \eqn{y} (given \eqn{\mu})
#'   for each generalized linear model family is given below:
#'   \itemize{
#'     \item{family: }{\eqn{Var(y)}}
#'     \item{poisson: }{\eqn{\mu \phi}}
#'     \item{nbinomial: }{\eqn{\mu + \mu^2 / \phi}}
#'     \item{binomial: }{\eqn{n \mu (1 - \mu) \phi}}
#'     \item{beta: }{\eqn{\mu (1 - \mu) / (1 + \phi)}}
#'     \item{Gamma: }{\eqn{\mu^2 / \phi}}
#'     \item{inverse.gaussian: }{\eqn{\mu^2 / \phi}}
#'   }
#'   The parameter \eqn{\phi} is a dispersion parameter that influences \eqn{Var(y)}.
#'   For the \code{poisson} and \code{binomial} families, \eqn{\phi} is always
#'   one. Note that this inverse Gaussian parameterization is different than a
#'   standard inverse Gaussian parameterization, which has variance \eqn{\mu^3 / \lambda}.
#'   Setting \eqn{\phi = \lambda / \mu} yields our parameterization, which is
#'   preferred for computational stability. Also note that the dispersion parameter
#'   is often defined in the literature as \eqn{V(\mu) \phi}, where \eqn{V(\mu)} is the variance
#'   function of the mean. We do not use this parameterization, which is important
#'   to recognize while interpreting dispersion parameter estimates using [spglm()] or [spgautor()].
#'   For more on generalized linear model constructions, see McCullagh and
#'   Nelder (1989).
#'
#' @return A list with two elements: \code{initial} and \code{is_known}.
#'   \code{initial} is a named numeric vector indicating the dispersion parameters
#'   with a specified initial and/or known value. \code{is_known} is a named
#'   numeric vector indicating whether the dispersion parameters in
#'   \code{initial} is known or not. The class of the list
#'   matches the value given to the \code{family} argument.
#'
#' @export
#'
#' @examples
#' # known dispersion value 1
#' dispersion_initial("nbinomial", dispersion = 1, known = "dispersion")
#' @references
#' McCullagh P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}. London: Chapman and Hall.
dispersion_initial <- function(family, dispersion, known) {

  # fix family
  if (is.symbol(substitute(family))) { # or is.language
    family <- deparse1(substitute(family))
  }

  # set defaults
  if (missing(dispersion)) dispersion <- NULL

  dispersion_checks(family, dispersion)

  dispersion_params_given <- c(dispersion = unname(dispersion))

  if (missing(known)) {
    is_known <- rep(FALSE, length(dispersion_params_given))
  } else {
    if (identical(known, "given")) {
      is_known <- rep(TRUE, length(dispersion_params_given))
    } else {
      is_known <- names(dispersion_params_given) %in% known
    }
  }
  names(is_known) <- names(dispersion_params_given)

  # error if NA and known
  dispersion_NA <- which(is.na(dispersion_params_given))
  if (any(is_known[dispersion_NA])) {
    stop("dispersion_initial values cannot be NA and known.", call. = FALSE)
  }

  new_dispersion_initial <- structure(list(initial = dispersion_params_given, is_known = is_known),
                                      class = family)
  new_dispersion_initial

}
