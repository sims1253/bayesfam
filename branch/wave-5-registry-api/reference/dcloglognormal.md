# Probability Density function of the Cloglognormal-Distribution

Probability Density function of the Cloglognormal-Distribution

## Usage

``` r
dcloglognormal(x, mu, sigma, log = FALSE)
```

## Arguments

- x:

  Value space of the function, x e (0, 1)

- mu:

  Median parameter, mu is already Cloglog-transformed, mu unbound

- sigma:

  Shape parameter, sigma \>= 0

- log:

  optional argument. If true, returns logarithmic probability. Default =
  FALSE

## Value

Normal Distribution density with Cloglog link function

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 1000)
plot(x, dcloglognormal(x, mu = 0.5, sigma = 2), type = "l")
```
