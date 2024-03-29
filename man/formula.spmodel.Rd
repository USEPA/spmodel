% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formula.R, R/formula_glm.R
\name{formula.spmodel}
\alias{formula.spmodel}
\alias{formula.splm}
\alias{formula.spautor}
\alias{formula.spglm}
\alias{formula.spgautor}
\title{Model formulae}
\usage{
\method{formula}{splm}(x, ...)

\method{formula}{spautor}(x, ...)

\method{formula}{spglm}(x, ...)

\method{formula}{spgautor}(x, ...)
}
\arguments{
\item{x}{A fitted model object from \code{\link[=splm]{splm()}}, \code{\link[=spautor]{spautor()}}, \code{\link[=spglm]{spglm()}}, or \code{\link[=spgautor]{spgautor()}}.}

\item{...}{Other arguments. Not used (needed for generic consistency).}
}
\value{
The formula used by a fitted model object.
}
\description{
Return formula used by a fitted model object.
}
\examples{
spmod <- splm(z ~ water + tarp,
  data = caribou,
  spcov_type = "exponential", xcoord = x, ycoord = y
)
formula(spmod)
}
