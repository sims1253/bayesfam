% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{inv_cauchit}
\alias{inv_cauchit}
\title{Cauchit response function}
\usage{
inv_cauchit(x)
}
\arguments{
\item{x}{value of x to be transformed, any real scalar or vector allowed}
}
\value{
inverse cauchit of x, result is e (0, 1)
}
\description{
Cauchit response function
}
\examples{
x <- seq(from = -10, to = 10, length.out = 100)
plot(x, inv_cauchit(x), type = "l")
}
