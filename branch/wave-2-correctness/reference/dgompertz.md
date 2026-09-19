# Probability density function for the Gompertz distribution, with Median parametrization.

Probability density function for the Gompertz distribution, with Median
parametrization.

## Usage

``` r
dgompertz(x, mu, beta, log = FALSE)
```

## Arguments

- x:

  Value

- mu:

  Median parameter

- beta:

  Scale parameter

- log:

  Optional argument. If TRUE, returns log(pdf).

## Value

f(x \| mu, eta)

## Details

PDF of Gompertz implementation, with constant b: \$\$b(\mu,\eta) := (1 /
\mu) \* log1p((-1 / \eta) \* log(0.5))\$\$ \$\$f(x) =
\eta\*b\*exp(\eta + bx - \eta \* e^{bx})\$\$

## Examples

``` r
x <- seq(from = 0.1, to = 5, length.out = 100)
plot(x, dgompertz(x, mu = 2, beta = 4), type = "l")
```
