% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sprnbinom.R
\name{sprnbinom}
\alias{sprnbinom}
\title{Simulate a spatial negative binomial random variable}
\usage{
sprnbinom(
  spcov_params,
  dispersion = 1,
  mean = 0,
  samples = 1,
  data,
  randcov_params,
  partition_factor,
  ...
)
}
\arguments{
\item{spcov_params}{An \code{\link[=spcov_params]{spcov_params()}} object.}

\item{dispersion}{The dispersion value.}

\item{mean}{A numeric vector representing the mean. \code{mean} must have length 1
(in which case it is recycled) or length equal
to the number of rows in \code{data}. The default is \code{0}.}

\item{samples}{The number of independent samples to generate. The default
is \code{1}.}

\item{data}{A data frame or \code{sf} object containing spatial information.}

\item{randcov_params}{A \code{\link[=randcov_params]{randcov_params()}} object.}

\item{partition_factor}{A formula indicating the partition factor.}

\item{...}{Additional arguments passed to \code{\link[=sprnorm]{sprnorm()}}.}
}
\value{
If \code{samples} is 1, a vector of random variables for each row of \code{data}
is returned. If \code{samples} is greater than one, a matrix of random variables
is returned, where the rows correspond to each row of \code{data} and the columns
correspond to independent samples.
}
\description{
Simulate a spatial negative binomial random variable with a specific
mean and covariance structure.
}
\details{
The values of \code{spcov_params}, \code{mean}, and \code{randcov_params}
are assumed to be on the link scale. They are used to simulate a latent normal (Gaussian)
response variable using \code{\link[=sprnorm]{sprnorm()}}. This latent variable is the
conditional mean used with \code{dispersion} to simulate a negative binomial random variable.
}
\examples{
spcov_params_val <- spcov_params("exponential", de = 0.2, ie = 0.1, range = 1)
sprnbinom(spcov_params_val, data = caribou, xcoord = x, ycoord = y)
sprnbinom(spcov_params_val, samples = 5, data = caribou, xcoord = x, ycoord = y)
}
