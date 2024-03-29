% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lomax.R
\name{lomax}
\alias{lomax}
\title{Lomax Stan-implementation in Mean parametrization.}
\usage{
lomax(link = "log", link_alpha = "log1p")
}
\arguments{
\item{link}{Link function for function}

\item{link_alpha}{Link function for eta argument}
}
\value{
brms Lomax distribution family
}
\description{
Lomax Stan-implementation in Mean parametrization.
}
\examples{
a <- rnorm(1000)
data <- list(a = a, y = rlomax(1000, exp(0.5 * a + 1), 2))
fit <- brms::brm(
  formula = y ~ 1 + a, data = data,
  family = lomax(), stanvars = lomax()$stanvars,
  refresh = 0
)
plot(fit)
}
