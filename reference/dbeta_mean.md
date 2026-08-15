# Mean parameterization the beta pdf.

Mean parameterization the beta pdf.

## Usage

``` r
dbeta_mean(x, mu, phi, log = FALSE)
```

## Arguments

- x:

  x-value, x e (0, 1)

- mu:

  Mean parameter, mu e (0, 1)

- phi:

  Precision parameter, phi \> 0

- log:

  Optional argument. If TRUE, returns log(pdf). Normally False.

## Value

PDF of custom Beta Distribution

## Details

The Beta Distribution has Density \$\$f(y \| \mu, \phi) =
\frac{\Gamma(\phi) x^{\mu\phi - 1} (1 -
x)^{(1-\mu)\phi}}{\Gamma(\mu\phi)\Gamma((1 - \mu)\phi)} \$\$

With parameterisation of the usual Beta-Distribution's shape parameters
a and b as: \$\$a := \mu\phi, b := (1 - \mu)\phi\$\$

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 1000)
plot(x, dbeta_mean(x, mu = 0.5, phi = 1), type = "l")
```
