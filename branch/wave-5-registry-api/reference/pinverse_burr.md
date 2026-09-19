# Cumulative distribution function for the Inverse Burr distribution, with Median parametrization.

Cumulative distribution function for the Inverse Burr distribution, with
Median parametrization.

## Usage

``` r
pinverse_burr(q, mu, tau, gamma)
```

## Source

Parameterization follows the actuar package
(<https://search.r-project.org/CRAN/refmans/actuar/html/InverseBurr.html>).

## Arguments

- q:

  Value to evaluate the CDF at, q \> 0

- mu:

  Median parameter of pdf, mu \> 0

- tau:

  Shape1 parameter of pdf, tau \> 0 (actuar: shape1)

- gamma:

  Shape2 parameter of pdf, gamma \> 0 (actuar: shape2)

## Value

CDF of the Inverse Burr distribution

## Details

Define the scale parameter sigma as \$\$\sigma(\mu, \tau, \gamma) := \mu
\cdot (2^{1/\tau} - 1)^{1/\gamma}\$\$

\$\$F(y \| \mu, \tau, \gamma) = (\frac{(y/\sigma)^{\gamma}}{1 +
(y/\sigma)^{\gamma}})^{\tau}\$\$

## Examples

``` r
x <- seq(from = 0.1, to = 10, length.out = 100)
plot(x, pinverse_burr(x, mu = 2, tau = 2, gamma = 3), type = "l")
```
