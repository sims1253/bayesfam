# Probability density function for the beta prime distribution (aka. inverse Beta)

Probability density function for the beta prime distribution (aka.
inverse Beta)

## Usage

``` r
dbetaprime(x, mu, phi, log = FALSE)
```

## Source

Bases on Bourguignon, M., Santos-Neto, M., & de Castro, M. (2018). A new
regression model for positive data (<https://arxiv.org/abs/1804.07734>)

## Arguments

- x:

  Value, x \> 0.

- mu:

  Mean, mu \> 0.

- phi:

  Precision, phi \> 0.

- log:

  Optional argument. If true, returns the log density.

## Value

Density of the pdf given x, mu and phi.

## Details

The beta prime distribution has density \$\$f(y) =
\frac{y^{(\mu(\Phi+1)-1)} (1+y)^{(-(\mu(\Phi+1)+\Phi+2))}}
{\Beta(\mu(1+\Phi), \Phi +2)}\$\$

With the usual beta prime parameters \$\$\beta = \Phi + 2, \alpha =
\mu(\Phi + 1)\$\$

## Examples

``` r
x <- seq(from = 0.1, to = 20, length.out = 1000)
plot(x, dbetaprime(x, mu = 4, phi = 2), type = "l")
```
