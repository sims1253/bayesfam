# Softplus density distribution in median parametrization.

Softplus density distribution in median parametrization.

## Usage

``` r
dsoftplusnormal(x, mu, sigma, log = FALSE)
```

## Arguments

- x:

  Value space of the distribution, x \> 0

- mu:

  Median parameter, mu is already log-transformed, mu unbound

- sigma:

  Sigma shape parameter, sigma \>= 0

- log:

  Bool argument, if true, returns the logarithmic density

## Value

Normal distribution density with logit link function

## Examples

``` r
x <- seq(from = 0.01, to = 10, length.out = 1000)
plot(x, dsoftplusnormal(x, mu = 1, sigma = 2), type = "l")
```
