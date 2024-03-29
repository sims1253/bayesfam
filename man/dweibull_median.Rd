% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weibull_median.R
\name{dweibull_median}
\alias{dweibull_median}
\title{Median parameterization of the Weibull pdf.}
\usage{
dweibull_median(x, mu, k, log = FALSE)
}
\arguments{
\item{x}{Value space, x > 0.}

\item{mu}{Median parameter, mu > 0.}

\item{k}{Shape parameter, k > 0.}

\item{log}{Optional argument. If TRUE, returns log(pdf). Normally False.}
}
\value{
f(x | mu, k)
}
\description{
Median parameterization of the Weibull pdf.
}
\details{
Define constant sigma as
\deqn{\sigma(\mu, k) := \mu / \Gamma(1 + 1 / k)}

The Weibull distribution density is defined as
\deqn{f(y) = \frac{k}{\sigma} * (\frac{x}{\sigma})^{\alpha - 1} * exp(-(\frac{x}{\sigma})^\alpha)}
}
\examples{
x <- seq(from = 0.01, to = 10, length.out = 1000)
plot(x, dweibull_median(x, mu = 2, k = 1), type = "l")
}
