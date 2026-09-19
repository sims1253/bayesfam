# Quantile function of Unit-Lindley distribution

Quantile function of Unit-Lindley distribution

## Usage

``` r
qunit_lindley(p, mu)
```

## Arguments

- p:

  vector of probabilities

- mu:

  Mean, mu e (0, 1)

## Value

q(p \| mu)

## Examples

``` r
p <- seq(from = 0.1, to = 0.9, length.out = 100)
plot(p, qunit_lindley(p, mu = 0.5))
```
