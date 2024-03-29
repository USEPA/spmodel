% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sprnorm.R
\name{sprnorm}
\alias{sprnorm}
\alias{sprnorm.exponential}
\alias{sprnorm.none}
\alias{sprnorm.car}
\title{Simulate a spatial normal (Gaussian) random variable}
\usage{
sprnorm(
  spcov_params,
  mean = 0,
  samples = 1,
  data,
  randcov_params,
  partition_factor,
  ...
)

\method{sprnorm}{exponential}(
  spcov_params,
  mean = 0,
  samples = 1,
  data,
  randcov_params,
  partition_factor,
  xcoord,
  ycoord,
  ...
)

\method{sprnorm}{none}(
  spcov_params,
  mean = 0,
  samples = 1,
  data,
  randcov_params,
  partition_factor,
  ...
)

\method{sprnorm}{car}(
  spcov_params,
  mean = 0,
  samples = 1,
  data,
  randcov_params,
  partition_factor,
  W,
  row_st = TRUE,
  M,
  ...
)
}
\arguments{
\item{spcov_params}{An \code{\link[=spcov_params]{spcov_params()}} object.}

\item{mean}{A numeric vector representing the mean. \code{mean} must have length 1
(in which case it is recycled) or length equal
to the number of rows in \code{data}. The default is \code{0}.}

\item{samples}{The number of independent samples to generate. The default
is \code{1}.}

\item{data}{A data frame or \code{sf} object containing spatial information.}

\item{randcov_params}{A \code{\link[=randcov_params]{randcov_params()}} object.}

\item{partition_factor}{A formula indicating the partition factor.}

\item{...}{Other arguments. Not used (needed for generic consistency).}

\item{xcoord}{Name of the column in \code{data} representing the x-coordinate.
Can be quoted or unquoted. Not required if \code{data} are an \code{sf}
object.}

\item{ycoord}{Name of the column in \code{data} representing the y-coordinate.
Can be quoted or unquoted. Not required if \code{data} are an \code{sf}
object.}

\item{W}{Weight matrix specifying the neighboring structure used for car and
sar models. Not required if \code{data} are an \code{sf}
polygon object and \code{W} should be calculated internally (using queen contiguity).}

\item{row_st}{A logical indicating whether row standardization be performed on
\code{W}. The default is \code{TRUE}.}

\item{M}{M matrix satisfying the car symmetry condition. The car
symmetry condition states that \eqn{(I - range * W)^{-1}M} is symmetric, where
\eqn{I} is an identity matrix, \eqn{range} is a constant that controls the
spatial dependence, \code{W} is the weights matrix,
and \eqn{^{-1}} represents the inverse operator.
\code{M} is required for car models
when \code{W} is provided and \code{row_st} is \code{FALSE}.  When \code{M},
is required, the default is the identity matrix.}
}
\value{
If \code{samples} is 1, a vector of random variables for each row of \code{data}
is returned. If \code{samples} is greater than one, a matrix of random variables
is returned, where the rows correspond to each row of \code{data} and the columns
correspond to independent samples.
}
\description{
Simulate a spatial normal (Gaussian) random variable with a specific
mean and covariance structure.
}
\details{
Random variables are simulated via the product of the covariance matrix's
square (Cholesky) root and independent standard normal random variables
with mean 0 and variance 1. Computing the square root is a significant
computational burden and likely unfeasible for sample sizes much past 10,000.
Because this square root only needs to be computed once, however, it is
nearly the sample computational cost to call \code{sprnorm()} for any value
of \code{samples}.

Only methods for the \code{exponential}, \code{none}, and \code{car}
covariance functions are documented here,
but methods exist for all other spatial covariance functions defined in
\code{\link[=spcov_initial]{spcov_initial()}}. Syntax for the \code{exponential} method is the same
as syntax for \code{spherical}, \code{gaussian}, \code{triangular},
\code{circular}, \code{cubic}, \code{pentaspherical}, \code{cosine}, \code{wave},
\code{jbessel}, \code{gravity}, \code{rquad}, \code{magnetic}, \code{matern},
\code{cauchy}, and \code{pexponential} methods. Syntax for
the \code{car} method is the same as syntax for the \code{sar} method. The
\code{extra} parameter for car and sar models is ignored when all observations have
neighbors.
}
\examples{
spcov_params_val <- spcov_params("exponential", de = 1, ie = 1, range = 1)
sprnorm(spcov_params_val, data = caribou, xcoord = x, ycoord = y)
sprnorm(spcov_params_val, mean = 1:30, samples = 5, data = caribou, xcoord = x, ycoord = y)
}
