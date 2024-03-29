% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unitlindley.R
\name{qunit_lindley}
\alias{qunit_lindley}
\title{Quantile function of Unit-Lindley distribution}
\usage{
qunit_lindley(p, mu)
}
\arguments{
\item{p}{vector of probabilities}

\item{mu}{Mean, mu e (0, 1)}
}
\value{
q(p | mu)
}
\description{
Quantile function of Unit-Lindley distribution
}
\examples{
p <- seq(from = 0.1, to = 0.9, length.out = 100)
plot(p, qunit_lindley(p, mu = 0.5))
}
