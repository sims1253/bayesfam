% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lomax.R
\name{qlomax}
\alias{qlomax}
\title{Quantile function for the Lomax distribution, with Mean parametrization.}
\usage{
qlomax(p, mu, alpha)
}
\arguments{
\item{p}{Quantile to be calculated}

\item{mu}{Median argument of Lomax}

\item{alpha}{Alpha argument of Gompertz}
}
\value{
Inverse of CDF, calculates a value, given a probability p
}
\description{
Quantile function for the Lomax distribution, with Mean parametrization.
}
\examples{
x <- seq(from = 0, to = 1, length.out = 100)
plot(x, qlomax(x, mu = 1, alpha = 2), type = "l")
}
