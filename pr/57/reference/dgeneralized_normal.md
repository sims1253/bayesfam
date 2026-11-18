# Probability density function for the generalized_normal distribution

Probability density function for the generalized_normal distribution

## Usage

``` r
dgeneralized_normal(x, mu, sigma, beta, log = FALSE)
```

## Arguments

- x:

  Value, unbound

- mu:

  Mean, unbound

- sigma:

  Scale, sigma \> 0

- beta:

  shape, beta \> 0

- log:

  Optional argument. If true, returns the log density.

## Value

Density of the pdf given x, mu and sigma

## Details

The generalized_normal distribution has density \$\$f(y \| \mu, \sigma,
\beta) = \frac{\beta}{2 \beta \Gamma(1/\beta)}exp(-\|z\|^\beta)\$\$

Where z is the linear transformation \$\$z(y, \mu, \sigma) = \frac{y -
\mu}{\sigma}\$\$

## Examples

``` r
x <- seq(from = -5, to = 10, length.out = 1000)
plot(x, dgeneralized_normal(x, mu = 1, sigma = 1, beta = 0.5), type = "l")
```
