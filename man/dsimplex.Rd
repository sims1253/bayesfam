% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplex.R
\name{dsimplex}
\alias{dsimplex}
\title{Simplex density function in mean parametrisation.}
\usage{
dsimplex(x, mu, sigma, log = FALSE)
}
\arguments{
\item{x}{value space, x e (0, 1)}

\item{mu}{Median parameter of pdf, mu e (0, 1)}

\item{sigma}{shape parameter, sigma unbound}

\item{log}{if true, returns log(pdf). Normally FALSE.}
}
\value{
f(x | mu, sigma)
}
\description{
Simplex density function in mean parametrisation.
}
\details{
\deqn{f(y) = (2 \pi \sigma^2(y(1-y))^3)^{-\frac{1}{2}} exp(-(\frac{y-\mu}{\mu(1-\mu)})^2 \frac{1}{2y(1-y)\sigma^2} )}
}
\examples{
x <- seq(from = 0.01, to = 0.99, length.out = 1000)
plot(x, dsimplex(x, mu = 0.7, sigma = 2), type = "l")
}
