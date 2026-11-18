# Density of the Unit Lindley likelihood

Density of the Unit Lindley likelihood

## Usage

``` r
dunit_lindley(x, mu, log = FALSE)
```

## Arguments

- x:

  vector of values, x e (0, 1)

- mu:

  Mean, mu e (0, 1)

- log:

  logical; if TRUE, log(pdf) is returned

## Value

f(x \| mu)

## Examples

``` r
x <- seq(from = 0.1, to = 0.9, length.out = 1000)
plot(x, dunit_lindley(x, 0.5))
```
