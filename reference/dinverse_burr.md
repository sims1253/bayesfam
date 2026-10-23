# Probability density function for the Inverse Burr distribution, with Median parametrization.

Probability density function for the Inverse Burr distribution, with
Median parametrization.

## Usage

``` r
dinverse_burr(x, mu, tau, gamma, log = FALSE)
```

## Source

Parameterization follows the actuar package (also known as the Dagum
distribution), with \\\tau\\ = shape1, \\\gamma\\ = shape2 and scale
\\\sigma = \mu (2^{1/\tau} - 1)^{1/\gamma}\\.

## Arguments

- x:

  Model space, defined for x \> 0

- mu:

  Median parameter of pdf, mu \> 0

- tau:

  Shape1 parameter of pdf, tau \> 0

- gamma:

  Shape2 parameter of pdf, gamma \> 0

- log:

  Optional log argument, if true, return log(pdf)

## Value

PDF of Inverse Burr Distribution

## Details

\$\$f(y \| \mu, \tau, \gamma) = \frac{\tau \gamma (y/\sigma)^{\gamma
\tau}}{y \[1 + (y/\sigma)^{\gamma}\]^{\tau + 1}}\$\$

## Examples

``` r
x <- seq(from = 0.1, to = 10, length.out = 100)
plot(x, dinverse_burr(x, mu = 2, tau = 2, gamma = 3), type = "l")
```
