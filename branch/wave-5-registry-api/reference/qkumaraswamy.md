# Quantile function of the Kumaraswamy distribution in Median parametrisation.

Quantile function of the Kumaraswamy distribution in Median
parametrisation.

## Usage

``` r
qkumaraswamy(u, mu = 0.5, p = 1)
```

## Arguments

- u:

  Quantile to be calculated, u e (0, 1)

- mu:

  Median parameter, mu e (0, 1)

- p:

  Phi shape parameter, p \> 0

## Value

q(u \| mu, p)

## Examples

``` r
u <- seq(from = 0.01, to = 0.09, length.out = 1000)
plot(u, qkumaraswamy(u, mu = 0.5, p = 2), type = "l")
```
