% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randcov_initial.R
\name{randcov_initial}
\alias{randcov_initial}
\title{Create a random effects covariance parameter initial object}
\usage{
randcov_initial(..., known)
}
\arguments{
\item{...}{Arguments to \code{randcov_params()}.}

\item{known}{A character vector indicating which random effect variances are to be
assumed known. The value \code{"given"} is shorthand for assuming all
random effect variances given to \code{randcov_initial()} are assumed known.}
}
\value{
A list with two elements: \code{initial} and \code{is_known}.
\code{initial} is a named numeric vector indicating the random effect variances
with specified initial and/or known values. \code{is_known} is a named
logical vector indicating whether the random effect variances in
\code{initial} are known or not.
}
\description{
Create a random effects (co)variance parameter initial object that specifies
initial and/or known values to use while estimating random effect variances
with modeling functions.
}
\details{
A random effect is specified as \eqn{Zu}, where \eqn{Z} is the random
effects design matrix and \code{u} is the random effect. The covariance of
\eqn{Zu} is \eqn{\sigma 2 ZZ^T}, where \eqn{\sigma 2} is the random effect
variance, and \eqn{Z^T} is the transpose of \eqn{Z}.
}
\examples{
randcov_initial(group = 1)
randcov_initial(group = 1, known = "group")
}
