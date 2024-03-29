% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-helper.R
\name{test_rng_quantiles}
\alias{test_rng_quantiles}
\title{Tests if an RNG can recover the true quantiles within a margin of error}
\usage{
test_rng_quantiles(
  rng_fun,
  quantile_fun,
  n,
  mu_list,
  aux_list = NA,
  aux2_list = NA,
  eps,
  quantiles,
  p_acceptable_failures,
  mu_link = identity,
  relative = FALSE
)
}
\arguments{
\item{rng_fun}{RNG function under test}

\item{quantile_fun}{Quantile function related to the rng under test.}

\item{n}{Sample size for the rng test.}

\item{mu_list}{Metric data used as RNG argument and to be compared to
(usually mean or median)}

\item{aux_list}{Auxiliary parameter value list.}

\item{aux2_list}{Auxiliary parameter value list for applicable distributions.}

\item{eps}{Acceptable difference of |mu - metric_mu(rng_fun)}

\item{quantiles}{Quantiles to test for recovery.}

\item{p_acceptable_failures}{Acceptable rate of failure, relative value of
difference bigger mu_eps}

\item{mu_link}{Default=identity, optional link-function argument, for example
useful in link-normal-distributions}

\item{relative}{True if the error should be relative to the mu_list}
}
\value{
Nothing actually, just wraps the test
}
\description{
Tests if an RNG can recover the true quantiles within a margin of error
}
\examples{
eps <- 0.001
mu_list <- seq(from = 1 + eps, to = 20, length.out = 10)
phi_list <- seq(from = 2 + eps, to = 20, length.out = 10)
# if working as expected, this test should not print any errors
bayesfam:::test_rng_quantiles(
  rng_fun = rbetaprime,
  quantile_fun = qbetaprime,
  n = 10000,
  mu_list = mu_list,
  aux_list = phi_list,
  eps = 0.1,
  quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
  p_acceptable_failures = 0.1,
  relative = TRUE
)
}
