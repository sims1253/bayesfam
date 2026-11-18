# Kumaraswamy density function in median parametrisation.

Kumaraswamy density function in median parametrisation.

## Usage

``` r
dkumaraswamy(x, mu, p, log = FALSE)
```

## Arguments

- x:

  vector of quantiles, x e (0, 1)

- mu:

  Median, mu e (0, 1)

- p:

  shape, p \> 0

- log:

  logical; if TRUE, log(pdf) is returned

## Value

f(x \| mu, p)

## Details

\$\$q(\mu, p) = -\frac{log(2)}{log(1-\mu^p)}\$\$

\$\$f(y \| \mu, p) = pqx^{p-1}(1-x^p)^{q-1}\$\$

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 1000)
plot(x, dkumaraswamy(x, mu = 0.5, p = 2), type = "l")
```
