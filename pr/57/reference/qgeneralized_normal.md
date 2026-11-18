# Quantile function of the generalized_normal distribution

Quantile function of the generalized_normal distribution

## Usage

``` r
qgeneralized_normal(p, mu, sigma, beta)
```

## Arguments

- p:

  quantile value, 0 \< p \< 1

- mu:

  Mean, unbound

- sigma:

  Scale, sigma \> 0

- beta:

  shape, beta \> 0

## Value

Quantiles of the generalized_normal distribution

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 100)
plot(x, qgeneralized_normal(x, mu = 2, sigma = 2, beta = 1), type = "l")
```
