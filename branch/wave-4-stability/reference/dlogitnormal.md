# Logitnormal density distribution in median parametrization.

Logitnormal density distribution in median parametrization.

## Usage

``` r
dlogitnormal(x, mu, sigma, log = FALSE)
```

## Arguments

- x:

  Value space of the distribution, x e (0, 1)

- mu:

  Median parameter, mu is already logit-transformed, mu unbound

- sigma:

  Sigma shape parameter, sigma \>= 0

- log:

  Bool argument, if true, returns the logarithmic density

## Value

Normal distribution density with logit link function

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 1000)
plot(x, dlogitnormal(x, mu = 0.5, sigma = 2), type = "l")
```
