% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cloglognormal.R
\name{rcloglognormal}
\alias{rcloglognormal}
\title{Cloglognormal RNG-function}
\usage{
rcloglognormal(n, mu = -0.36, sigma = 0.75)
}
\arguments{
\item{n}{Number of draws}

\item{mu}{Median parameter, mu unbound, mu already cloglog transformed}

\item{sigma}{Shape parameter}
}
\value{
n cloglog-normally distributed samples
}
\description{
Cloglognormal RNG-function
}
\examples{
hist(rcloglognormal(100, 0.5, 2))
}
