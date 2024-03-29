% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/softplusnormal.R
\name{dsoftplusnormal}
\alias{dsoftplusnormal}
\title{Softplus density distribution in median parametrization.}
\usage{
dsoftplusnormal(x, mu, sigma, log = FALSE)
}
\arguments{
\item{x}{Value space of the distribution, x > 0}

\item{mu}{Median parameter, mu is already log-transformed, mu unbound}

\item{sigma}{Sigma shape parameter, sigma >= 0}

\item{log}{Bool argument, if true, returns the logarithmic density}
}
\value{
Normal distribution density with logit link function
}
\description{
Softplus density distribution in median parametrization.
}
\examples{
x <- seq(from = 0.01, to = 10, length.out = 1000)
plot(x, dsoftplusnormal(x, mu = 1, sigma = 2), type = "l")
}
