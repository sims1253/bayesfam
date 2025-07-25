% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-helper.R
\name{test_brms_quantile}
\alias{test_brms_quantile}
\title{Check, that data of the posterior is close enough to the reference data.}
\usage{
test_brms_quantile(posterior_data, arg_name, reference, thresh, debug = FALSE)
}
\arguments{
\item{posterior_data}{Data fitted and drawn, a brms object, which is a list in R terms}

\item{arg_name}{Name of the argument variable to check, as single string}

\item{reference}{Reference value to check against, single real scalar}

\item{thresh}{real scalar or 2-length vector of quantile bounds.
For scalar constructs bound as \verb{[thresh, 1-thresh]}
thresh has to be inside the Unit-Interval.}

\item{debug}{True for verbose output of test results.}
}
\value{
Single boolean success, fail or error
}
\description{
Check, that data of the posterior is close enough to the reference data.
}
\examples{
fit <- bayesfam:::construct_brms(
  n_data_sampels = 1000,
  intercept = 5.0,
  aux_par = 2.0,
  rng_link = identity,
  family = betaprime,
  rng = rbetaprime
)
result <- bayesfam:::test_brms_quantile(
  posterior_data = fit, arg_name = "phi", 2.0, 0.025
)
plot(fit)
# beta_prime uses log-link for Intercept
}
