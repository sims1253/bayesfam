# Mean parameterization of the beta quantile function

Mean parameterization of the beta quantile function

## Usage

``` r
qbeta_mean(p, mu, phi)
```

## Arguments

- p:

  Probabilities, for which to calculate the quantiles

- mu:

  Mean, mu e (0, 1).

- phi:

  Precision parameter, phi \> 0.

## Value

Quantiles of the beta distribution

## Examples

``` r
x <- seq(from = 0, to = 1, length.out = 100)
plot(x, qbeta_mean(x, mu = 0.5, phi = 2), type = "l")
```
