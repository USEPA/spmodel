% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{moose}
\alias{moose}
\title{Moose counts and presence in Alaska, USA}
\format{
An \code{sf} object with 218 rows and 5 columns.
\itemize{
\item elev: The elevation.
\item strat: A factor representing strata (used for sampling). Can take values \code{L} and \code{M}.
\item count: The count (number) of moose observed.
\item presence: A binary factor representing whether no moose were observed (value \code{0}) or at least one moose was observed
(value \code{1}).
\item geometry: \code{POINT} geometry representing coordinates in an Alaska
Albers projection (EPSG: 3338). Distances between points are in meters.
}
}
\source{
Alaska Department of Fish and Game, Division of Wildlife Conservation has released
this data set under the CC0 license.
}
\usage{
moose
}
\description{
Moose counts and presence in Alaska, USA.
}
\keyword{datasets}
