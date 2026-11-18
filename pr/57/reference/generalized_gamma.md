# Generalized gamma distribution

Generalized gamma distribution

## Usage

``` r
dgeneralized_gamma(x, mu = 0, sigma = 1, Q, log = FALSE)

rgeneralized_gamma(n, mu = 0, sigma = 1, Q)

generalized_gamma(link = "log", link_sigma = "log", link_Q = "log")
```

## Source

Bases on flexsurv (<https://github.com/chjackson/flexsurv/tree/master>)
by Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>. Inspired by a
blog post by Demetri Pananos
(<https://dpananos.github.io/posts/2023-12-02-gen-gamma/>) and code by
Krzysztof Sakrejda (<https://github.com/sakrejda/tooling>).

## Arguments

- x:

  Value, x \> 0.

- mu:

  Vector of “location” parameters.

- sigma:

  Vector of
  `scale'' parameters. Note the inconsistent meanings of the term `scale” -
  this parameter is analogous to the (log-scale) standard deviation of
  the log-normal distribution,
  `sdlog'' in [dlnorm()], rather than the `scale” parameter of the gamma
  distribution [`dgamma()`](https://rdrr.io/r/stats/GammaDist.html).
  Constrained to be positive.

- Q:

  Vector of shape parameters.

- log:

  logical; if TRUE the log-pdf is returned

- n:

  number of observations.

- link:

  Link function for mu

- link_sigma:

  Link function for sigma

- link_Q:

  Link function for Q

## Value

`dgeneralized_gamma` gives the density, `rgeneralized_gamma` generates
random deviates,

## References

Prentice, R. L. (1974). A log gamma model and its maximum likelihood
estimation. Biometrika 61(3):539-544.

Farewell, V. T. and Prentice, R. L. (1977). A study of distributional
shape in life testing. Technometrics 19(1):69-75.

Lawless, J. F. (1980). Inference in the generalized gamma and log gamma
distributions. Technometrics 22(3):409-419.

Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007). Parametric
survival analysis and taxonomy of hazard functions for the generalized
gamma distribution. Statistics in Medicine 26:4252-4374

Stacy, E. W. (1962). A generalization of the gamma distribution. Annals
of Mathematical Statistics 33:1187-92
