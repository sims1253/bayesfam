# Quantile function of the logistic distribution

Quantile function of the logistic distribution

## Usage

``` r
qlogistic(p, mu, sigma)
```

## Arguments

- p:

  quantile value, 0 \< p \< 1

- mu:

  Mean, unbound

- sigma:

  Scale, sigma \> 0

## Value

Quantiles of the logistic distribution

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 100)
plot(x, qlogistic(x, mu = 2, sigma = 2), type = "l")
```
