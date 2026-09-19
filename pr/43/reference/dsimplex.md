# Simplex density function in mean parametrisation.

Simplex density function in mean parametrisation.

## Usage

``` r
dsimplex(x, mu, sigma, log = FALSE)
```

## Arguments

- x:

  value space, x e (0, 1)

- mu:

  Median parameter of pdf, mu e (0, 1)

- sigma:

  shape parameter, sigma unbound

- log:

  if true, returns log(pdf). Normally FALSE.

## Value

f(x \| mu, sigma)

## Details

\$\$f(y) = (2 \pi \sigma^2(y(1-y))^3)^{-\frac{1}{2}}
exp(-(\frac{y-\mu}{\mu(1-\mu)})^2 \frac{1}{2y(1-y)\sigma^2} )\$\$

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 1000)
plot(x, dsimplex(x, mu = 0.7, sigma = 2), type = "l")
```
