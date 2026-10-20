# Quantile function of the gumbel distribution

Quantile function of the gumbel distribution

## Usage

``` r
qgumbel_mean(p, mu, sigma)
```

## Arguments

- p:

  quantile value, 0 \< p \< 1

- mu:

  Mean, unbound

- sigma:

  Scale, sigma \> 0

## Value

Quantiles of the gumbel distribution

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 100)
plot(x, qgumbel_mean(x, mu = 2, sigma = 2), type = "l")
```
