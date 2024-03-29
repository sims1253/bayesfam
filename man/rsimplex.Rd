% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplex.R
\name{rsimplex}
\alias{rsimplex}
\title{Simplex RNG function in Median parametrization.}
\usage{
rsimplex(n, mu = 0.5, sigma = 1)
}
\arguments{
\item{n}{Number of samples to draw, as a natural number scalar.}

\item{mu}{Mean parameter, mu e (0, 1)}

\item{sigma}{shape parameter, Sigma unbound}
}
\value{
n samples in Simplex distribution.
}
\description{
Based on code from simplexreg
Peng Zhang, Zhenguo Qiu, Chengchun Shi (2016). simplexreg: An R
Package for Regression Analysis of Proportional Data Using the
Simplex Distribution. Journal of Statistical Software, 71(11), 1-21.
doi:10.18637/jss.v071.i11
}
\examples{
hist(rsimplex(10000, mu = 0.7, sigma = 2))
}
