# Quantile function of the beta prime distribution

Quantile function of the beta prime distribution

## Usage

``` r
qbetaprime(p, mu, phi)
```

## Arguments

- p:

  Probabilities, for which to calculate the quantiles

- mu:

  Mean, mu \> 0.

- phi:

  Precision, phi \> 0.

## Value

Quantiles of the beta prime distribution

## Examples

``` r
x <- seq(from = 0, to = 1, length.out = 100)
plot(x, qbetaprime(x, mu = 1, phi = 2), type = "l")
```
