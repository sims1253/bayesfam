% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cloglognormal.R
\name{cloglognormal}
\alias{cloglognormal}
\title{Custom brms family Cloglog-Normal in median parametrization.}
\usage{
cloglognormal(link = "identity", link_sigma = "log")
}
\arguments{
\item{link}{Link function argument (as string) for Median argument. Left as identity!}

\item{link_sigma}{Link function argument (as string) for Shape argument}
}
\value{
Cloglog brms model-object
}
\description{
Custom brms family Cloglog-Normal in median parametrization.
}
\examples{
data <- rcloglognormal(1000, 0.5, 2)
# cloglognormal does not like values to close to the boundary
data <- limit_data(data, c(1e-12, 1 - 1e-12))
fit <- brms::brm(
  formula = y ~ 1, data = list(y = data),
  family = cloglognormal(), stanvars = cloglognormal()$stanvars,
  refresh = 0
)
plot(fit)
}
