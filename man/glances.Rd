% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glances.R, R/glances_glm.R
\name{glances}
\alias{glances}
\alias{glances.splm}
\alias{glances.spautor}
\alias{glances.splm_list}
\alias{glances.spautor_list}
\alias{glances.spglm}
\alias{glances.spgautor}
\alias{glances.spglm_list}
\alias{glances.spgautor_list}
\title{Glance at many fitted model objects}
\usage{
glances(object, ...)

\method{glances}{splm}(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE)

\method{glances}{spautor}(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE)

\method{glances}{splm_list}(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE)

\method{glances}{spautor_list}(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE)

\method{glances}{spglm}(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE)

\method{glances}{spgautor}(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE)

\method{glances}{spglm_list}(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE)

\method{glances}{spgautor_list}(object, ..., sort_by = "AICc", decreasing = FALSE, warning = TRUE)
}
\arguments{
\item{object}{A fitted model object from \code{\link[=splm]{splm()}}, \code{\link[=spautor]{spautor()}}, \code{\link[=spglm]{spglm()}}, or \code{\link[=spgautor]{spgautor()}}.}

\item{...}{Additional fitted model objects. Ignored
if \code{object} has class \code{splm_list}, \code{spautor_list}, \code{spglm_list}, or \code{spgautor_list}.}

\item{sort_by}{Sort by a \code{glance} statistic (i.e., the name of a column
output from \code{glance()} or the order of model input (\code{sort_by = "order"}).
The default is \code{"AICc"}.}

\item{decreasing}{Whether \code{sort_by} should sort by decreasing order? The default is \code{FALSE}.}

\item{warning}{Whether a warning is displayed when model comparisons violate certain rules.
The default is \code{TRUE}.}
}
\value{
A tibble where each row represents the output of \code{glance()} for
each fitted model object.
}
\description{
\code{glances()} repeatedly calls \code{glance()} on several
fitted model objects and binds the output together, sorted by a column of interest.
}
\examples{
lmod <- splm(z ~ water + tarp,
  data = caribou,
  spcov_type = "none"
)
spmod <- splm(z ~ water + tarp,
  data = caribou,
  spcov_type = "exponential", xcoord = x, ycoord = y
)
glances(lmod, spmod)
glances(lmod, spmod, sort_by = "logLik", decreasing = TRUE)
}
