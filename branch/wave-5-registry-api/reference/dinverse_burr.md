# Probability density function for the Inverse Burr distribution, with Median parametrization.

Probability density function for the Inverse Burr distribution, with
Median parametrization.

## Usage

``` r
dinverse_burr(x, mu, tau, gamma, log = FALSE)
```

## Source

Parameterization follows the actuar package
(<https://search.r-project.org/CRAN/refmans/actuar/html/InverseBurr.html>).
Also known as the Dagum distribution.

## Arguments

- x:

  Model space, defined for x \> 0

- mu:

  Median parameter of pdf, mu \> 0

- tau:

  Shape1 parameter of pdf, tau \> 0 (actuar: shape1)

- gamma:

  Shape2 parameter of pdf, gamma \> 0 (actuar: shape2)

- log:

  Optional log argument, if true, return log(pdf)

## Value

PDF of the Inverse Burr distribution

## Details

Define the scale parameter sigma as \$\$\sigma(\mu, \tau, \gamma) := \mu
\cdot (2^{1/\tau} - 1)^{1/\gamma}\$\$

\$\$f(y \| \mu, \tau, \gamma) = \frac{\tau \gamma (y/\sigma)^{\gamma
\tau}}{y \[1 + (y/\sigma)^{\gamma}\]^{\tau + 1}}\$\$

## Examples

``` r
x <- seq(from = 0.1, to = 10, length.out = 100)
plot(x, dinverse_burr(x, mu = 2, tau = 2, gamma = 3), type = "l")
```
