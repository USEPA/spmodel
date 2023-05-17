#' Create a dispersion parameter object
#'
#' @description Create a dispersion parameter object for use with other
#'   functions.
#'
#' @param family The generalized linear model family describing the distribution
#'   of the response variable to be used. \code{"poisson"}, \code{"nbinomial"}, \code{"binomial"},
#'   \code{"beta"}, \code{"Gamma"}, and \code{"inverse.gaussian"}.
#' @param dispersion The value of the dispersion parameter.
#'
#' @details The variance function of an individual \eqn{y} (given \eqn{\mu})
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
#' @return A named numeric vector with class \code{family} containing the dispersion.
#'
#' @export
#'
#' @examples
#' dispersion_params("beta", dispersion = 1)#'
#' @references
#' McCullagh P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}. London: Chapman and Hall.
dispersion_params <- function(family, dispersion) {

  # fix family
  if (is.symbol(substitute(family))) { # or is.language
    family <- deparse1(substitute(family))
  }

  dispersion_checks(family, dispersion)
  # running after dispersion check
  if (! family %in% c("poisson", "nbinomial", "binomial", "beta", "Gamma", "inverse.gaussian")) {
    stop(paste(family, " is not a valid glm family for this function.", sep = ""), call. = FALSE)
  }

  object <- c(dispersion = unname(dispersion))
  new_object <- structure(object, class = family)
  new_object
}
