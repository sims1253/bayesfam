% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shifted_inv_gaussian.R
\name{rshifted_inv_gaussian}
\alias{rshifted_inv_gaussian}
\title{Shifted inverse Gaussian RNG function}
\usage{
rshifted_inv_gaussian(n, mu = 1, shape = 1, shift = 1)
}
\arguments{
\item{n}{number of observations}

\item{mu}{Mean, mu > 0}

\item{shape}{Shape, shape unbound}

\item{shift}{Shift, shift >= 0}
}
\value{
n samples drawn from the shifted inverse Gaussian distribution
}
\description{
Shifted inverse Gaussian RNG function
}
\examples{
hist(rshifted_inv_gaussian(100, 1, 1, 1))
}
