% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lognormal_natural.R
\name{rlognormal_natural}
\alias{rlognormal_natural}
\title{Lognormal Natural RNG function}
\usage{
rlognormal_natural(n, mu = 1, sigma = 1)
}
\arguments{
\item{n}{number of observations}

\item{mu}{mean, mu > 0}

\item{sigma}{sigma, sigma > 0}
}
\value{
n samples drawn from the Lognormal natural distribution
}
\description{
Lognormal Natural RNG function
}
\examples{
hist(rlognormal_natural(100, 1, 2))
}
