# Cauchitnormal Density distribution.

Cauchitnormal Density distribution.

## Usage

``` r
dcauchitnormal(x, mu, sigma, log = FALSE)
```

## Arguments

- x:

  Value space of the distribution, x e (0, 1)

- mu:

  Median parameter, mu is already Cauchit-transformed, mu unbound

- sigma:

  Sigma shape parameter, sigma \>= 0

- log:

  Bool argument, if true, returns the logarithmic density

## Value

Normal Distribution Density with Cauchit link function

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 1000)
plot(x, dcauchitnormal(x, mu = 0.5, sigma = 2), type = "l")
```
