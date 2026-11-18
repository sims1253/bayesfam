# Probability density function for the gumbel distribution

Probability density function for the gumbel distribution

## Usage

``` r
dgumbel_mean(x, mu, sigma, log = FALSE)
```

## Arguments

- x:

  Value, unbound

- mu:

  Mean, unbound

- sigma:

  Scale, sigma \> 0

- log:

  Optional argument. If true, returns the log density.

## Value

Density of the pdf given x, mu and sigma

## Details

The gumbel distribution has density \$\$f(y \| \mu, \sigma) =
\frac{1}{\sigma} exp(-(z + e^{-z}))\$\$

Where z is the linear transformation \$\$z(y, \mu, \sigma) = \frac{y -
\mu}{\sigma} + \gamma\$\$

Where Gamma refers to the Euler-Mascheroni constant

## Examples

``` r
x <- seq(from = -5, to = 10, length.out = 1000)
plot(x, dgumbel_mean(x, mu = 2, sigma = 2), type = "l")
```
